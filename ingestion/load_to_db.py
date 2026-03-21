"""
load_to_db.py
-------------

Calls transform.py, then loads the two clean DataFrames into PostgreSQL.
All three scripts (fetch → transform → load) can be run in sequence by
executing this file alone in the Docker container.

Load strategy:
    isolates  — upsert on biosample_acc   (idempotent re-runs)
    genes     — delete-then-insert for each batch of isolates
                (avoids duplicates; always reflects the latest NCBI data)
    summary   — TRUNCATE + repopulate from a GROUP BY across base tables
                (kept small for fast dashboard queries)
    pipeline_runs — insert on start, update on finish (full audit trail)

Environment variables:
    DB_HOST, DB_PORT, DB_NAME, DB_USER, DB_PASSWORD
    GCS_BUCKET, GCS_PREFIX, ORGANISM
    (all passed through to transform.py)
"""

import os
import logging
from datetime import datetime, UTC

import pandas as pd
import psycopg2
from psycopg2.extras import execute_values

from transform import main as run_transform

log = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

DB_CONFIG = {
    "host":     os.getenv("DB_HOST",     "localhost"),
    "port":     int(os.getenv("DB_PORT", "5432")),
    "dbname":   os.getenv("DB_NAME",     "amr_db"),
    "user":     os.getenv("DB_USER",     "amr_user"),
    "password": os.getenv("DB_PASSWORD", "amr_password"),
}

ORGANISM = os.getenv("ORGANISM", "Pseudomonas_aeruginosa")


# ── Connection ─────────────────────────────────────────────────────────────

def get_conn():
    return psycopg2.connect(**DB_CONFIG)


# ── Pipeline run tracking ──────────────────────────────────────────────────

def start_run(conn, organism: str, n_fetched: int) -> int:
    with conn.cursor() as cur:
        cur.execute(
            """
            INSERT INTO pipeline_runs (organism, records_fetched, status)
            VALUES (%s, %s, 'running') RETURNING run_id
            """,
            (organism, n_fetched),
        )
        run_id = cur.fetchone()[0]
    conn.commit()
    log.info("Pipeline run started → run_id=%d", run_id)
    return run_id


def finish_run(conn, run_id: int, n_loaded: int, status="success", notes=""):
    with conn.cursor() as cur:
        cur.execute(
            """
            UPDATE pipeline_runs
            SET records_loaded=%s, status=%s, notes=%s
            WHERE run_id=%s
            """,
            (n_loaded, status, notes, run_id),
        )
    conn.commit()
    log.info("Pipeline run %d → %s (%d records loaded)", run_id, status, n_loaded)

# ── Helper: safe row conversion ────────────────────────────────────────────

def safe_row_tuple(row, cols):
    return tuple(
        None if pd.isna(row[c]) else row[c]
        for c in cols
    )

# ── Upsert: isolates ───────────────────────────────────────────────────────

def upsert_isolates(conn, df, run_id: int) -> int:
    cols = [
        "biosample_acc", "organism", "strain",
        "collection_date", "collection_year",
        "geo_loc_name", "country",
        "isolation_source", "host", "serovar",
        "amr_genotypes", "run_id",
    ]
    for c in cols:
        if c not in df.columns:
            df[c] = None
    df["run_id"] = run_id

    rows = [safe_row_tuple(row, cols) for _, row in df.iterrows()]

    sql = f"""
        INSERT INTO isolate_metadata ({", ".join(cols)})
        VALUES %s
        ON CONFLICT (biosample_acc) DO UPDATE SET
            organism         = EXCLUDED.organism,
            collection_date  = EXCLUDED.collection_date,
            collection_year  = EXCLUDED.collection_year,
            geo_loc_name     = EXCLUDED.geo_loc_name,
            country          = EXCLUDED.country,
            amr_genotypes    = EXCLUDED.amr_genotypes,
            run_id           = EXCLUDED.run_id
    """
    with conn.cursor() as cur:
        execute_values(cur, sql, rows, page_size=500)
    conn.commit()
    log.info("Upserted %d isolate rows", len(rows))
    return len(rows)


# ── Replace: genes ─────────────────────────────────────────────────────────

def load_genes(conn, df) -> int:
    if df.empty:
        log.info("No gene rows to load.")
        return 0

    biosample_ids = list(df["biosample_acc"].unique())

    cols = [
        "biosample_acc", "gene_symbol", "element_type",
        "element_subtype", "drug_class", "subclass",
        "method", "scope", "pct_coverage", "pct_identity", "ref_seq_name",
    ]
    for c in cols:
        if c not in df.columns:
            df[c] = None

    rows = [safe_row_tuple(row, cols) for _, row in df.iterrows()]

    with conn.cursor() as cur:
        # Remove stale rows for this batch of isolates
        cur.execute(
            "DELETE FROM amr_genes WHERE biosample_acc = ANY(%s)",
            (biosample_ids,),
        )
        deleted = cur.rowcount

        execute_values(
            cur,
            f"INSERT INTO amr_genes ({', '.join(cols)}) VALUES %s",
            rows,
            page_size=1000,
        )
    conn.commit()
    log.info("Deleted %d stale gene rows; inserted %d fresh rows", deleted, len(rows))
    return len(rows)


# ── Refresh summary ────────────────────────────────────────────────────────

def refresh_summary(conn):
    """
    Repopulate amr_summary by aggregating across base tables.
    This table is the primary source for dashboard charts —
    keeping it pre-aggregated means dashboard queries stay fast
    even as the base tables grow to millions of rows.
    """
    sql = """
        BEGIN;
        TRUNCATE TABLE amr_summary;

        INSERT INTO amr_summary
            (organism, country, collection_year,
             drug_class, subclass, gene_symbol,
             isolate_count, gene_count)
        SELECT
            im.organism,
            im.country,
            im.collection_year,
            ag.drug_class,
            ag.subclass,
            ag.gene_symbol,
            COUNT(DISTINCT im.biosample_acc) AS isolate_count,
            COUNT(ag.id)                     AS gene_count
        FROM amr_genes ag
        JOIN isolate_metadata  im ON ag.biosample_acc = im.biosample_acc
        GROUP BY
            im.organism, im.country, im.collection_year,
            ag.drug_class, ag.subclass, ag.gene_symbol;

        COMMIT;
    """
    with conn.cursor() as cur:
        cur.execute(sql)
    log.info("amr_summary refreshed.")


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    log.info("=== AMR Load Pipeline starting ===")

    # Step 1 — Transform
    log.info("Running transform …")
    isolates_df, genes_df = run_transform()

    # Step 2 — Connect
    conn = get_conn()
    log.info("Connected to PostgreSQL at %s/%s", DB_CONFIG["host"], DB_CONFIG["dbname"])

    # Step 3 — Open run record
    run_id = start_run(conn, ORGANISM.replace("_", " "), len(isolates_df))

    try:
        n_iso   = upsert_isolates(conn, isolates_df, run_id)
        n_genes = load_genes(conn, genes_df)
        refresh_summary(conn)
        finish_run(conn, run_id, n_iso)

    except Exception as exc:
        conn.rollback() 
        finish_run(conn, run_id, 0, status="failed", notes=str(exc))
        conn.close()
        raise

    conn.close()
    log.info("=== AMR Load Pipeline complete ===")
    return {"run_id": run_id, "isolates_loaded": n_iso, "genes_loaded": n_genes}


if __name__ == "__main__":
    result = main()
    print(result)


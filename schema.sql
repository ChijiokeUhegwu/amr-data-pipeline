-- ============================================================
-- AMR Surveillance Pipeline — PostgreSQL Schema
-- Organism: Pseudomonas aeruginosa (configurable)
--
-- Table layout:
--   pipeline_runs     → audit log of every Kestra execution
--   isolate_metadata  → one row per isolate (geography, host, date)
--   amr_genes         → one row per detected AMR gene per isolate
--   amr_summary       → pre-aggregated for fast Streamlit queries
-- ============================================================

-- ── Run audit log ─────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS pipeline_runs (
    run_id          SERIAL PRIMARY KEY,
    run_date        TIMESTAMPTZ DEFAULT NOW(),
    organism        TEXT,
    records_fetched INTEGER,
    records_loaded  INTEGER,
    status          TEXT CHECK (status IN ('running', 'success', 'failed')),
    notes           TEXT
);

-- ── Isolate-level metadata ─────────────────────────────────────────────────
-- One row per NCBI biosample / isolate.
-- biosample_acc is the primary key and the foreign key used in amr_genes.
CREATE TABLE IF NOT EXISTS isolate_metadata (
    biosample_acc    TEXT PRIMARY KEY,
    organism         TEXT,
    strain           TEXT,
    collection_date  TEXT,           -- raw NCBI string, e.g. "2019-03" or "2019"
    collection_year  INTEGER,        -- extracted 4-digit year for easy aggregation
    geo_loc_name     TEXT,           -- raw NCBI field, e.g. "Nigeria: Lagos"
    country          TEXT,           -- parsed from geo_loc_name
    isolation_source TEXT,
    host             TEXT,
    serovar          TEXT,
    amr_genotypes    TEXT,           -- NCBI-provided text summary (informational)
    run_id           INTEGER REFERENCES pipeline_runs(run_id),
    inserted_at      TIMESTAMPTZ DEFAULT NOW()
);

-- ── AMR gene hits ──────────────────────────────────────────────────────────
-- One row per AMRFinderPlus hit per isolate.
-- Columns map directly to the NCBI AMR TSV after column name cleaning.
CREATE TABLE IF NOT EXISTS amr_genes (
    id               SERIAL PRIMARY KEY,
    biosample_acc    TEXT REFERENCES isolate_metadata(biosample_acc),

    gene_symbol      TEXT,           -- e.g. "blaOXA-50", "mexB"
    element_type     TEXT,           -- "AMR" (filtered), "STRESS", "VIRULENCE"
    element_subtype  TEXT,           -- gene family e.g. "BETA-LACTAM"
    drug_class       TEXT,           -- e.g. "BETA-LACTAM"
    subclass         TEXT,           -- e.g. "CARBAPENEM", "CEPHALOSPORIN"
    method           TEXT,           -- detection method: BLASTP, HMM, EXACTX, etc.
    scope            TEXT,           -- "core" or "plus"
    pct_coverage     NUMERIC(5,2),   -- % coverage of reference sequence
    pct_identity     NUMERIC(5,2),   -- % identity to reference sequence
    ref_seq_name     TEXT,           -- name of closest reference sequence

    inserted_at      TIMESTAMPTZ DEFAULT NOW()
);

-- ── Pre-aggregated summary ─────────────────────────────────────────────────
-- Repopulated on every pipeline run by load_to_db.py → refresh_summary().
-- Streamlit queries this table instead of aggregating the base tables live,
-- keeping dashboard response times fast even with large datasets.
CREATE TABLE IF NOT EXISTS amr_summary (
    id              SERIAL PRIMARY KEY,
    organism        TEXT,
    country         TEXT,
    collection_year INTEGER,
    drug_class      TEXT,
    subclass        TEXT,
    gene_symbol     TEXT,
    isolate_count   INTEGER,
    gene_count      INTEGER,
    last_updated    TIMESTAMPTZ DEFAULT NOW()
);

-- ── Indexes ────────────────────────────────────────────────────────────────
CREATE INDEX IF NOT EXISTS idx_isolate_country  ON isolate_metadata(country);
CREATE INDEX IF NOT EXISTS idx_isolate_year     ON isolate_metadata(collection_year);
CREATE INDEX IF NOT EXISTS idx_isolate_organism ON isolate_metadata(organism);

CREATE INDEX IF NOT EXISTS idx_gene_biosample   ON amr_genes(biosample_acc);
CREATE INDEX IF NOT EXISTS idx_gene_symbol      ON amr_genes(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_gene_class       ON amr_genes(drug_class);

CREATE INDEX IF NOT EXISTS idx_summary_country  ON amr_summary(country);
CREATE INDEX IF NOT EXISTS idx_summary_gene     ON amr_summary(gene_symbol);

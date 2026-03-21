"""
transform.py
------------

Reads the raw PDG*.amr.metadata.tsv from GCS (or a local file for dev)
and produces two clean DataFrames:

    isolates_df  -- one row per isolate (geography, host, date)
    genes_df     -- one row per AMR gene per isolate

The amr.metadata.tsv is a per-isolate summary file that already contains
both geographic metadata AND amr_genotypes (a comma-separated gene list).
Example amr_genotypes value:  "blaOXA-50,mexB,mexD,oprN,fosA"

This script:
  1. Selects + renames the geographic/clinical columns  -> isolates_df
  2. Explodes the amr_genotypes string into individual gene rows
  3. Classifies each gene into a drug_class / subclass   -> genes_df

Environment variables:
    GCS_BUCKET  -- bucket name
    GCS_PREFIX  -- e.g. raw/Pseudomonas_aeruginosa/20240601_120000
                   (auto-detected from latest run if not set)
    ORGANISM    -- used for auto-detection fallback
    INPUT_PATH  -- local TSV path for offline testing (skips GCS entirely)
"""

import io
import os
import re
import logging

import pandas as pd
from google.cloud import storage


log = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

GCS_BUCKET = os.getenv("GCS_BUCKET", "amr-data-pipeline-data-lake")
GCS_PREFIX = os.getenv("GCS_PREFIX", "")
ORGANISM   = os.getenv("ORGANISM",   "Pseudomonas_aeruginosa")
INPUT_PATH = os.getenv("INPUT_PATH", "")


# -- Load raw file ----------------------------------------------------------

def _latest_gcs_prefix() -> str:
    client  = storage.Client()
    blobs   = list(client.list_blobs(GCS_BUCKET, prefix=f"raw/{ORGANISM}/"))
    if not blobs:
        raise FileNotFoundError(f"No blobs under gs://{GCS_BUCKET}/raw/{ORGANISM}/")
    prefixes = sorted({b.name.rsplit("/", 1)[0] for b in blobs}, reverse=True)
    log.info("Auto-detected GCS prefix: %s", prefixes[0])
    return prefixes[0]


def load_raw() -> bytes:
    if INPUT_PATH:
        log.info("Dev mode: loading local file %s", INPUT_PATH)
        with open(INPUT_PATH, "rb") as f:
            return f.read()
    prefix = GCS_PREFIX or _latest_gcs_prefix()
    client = storage.Client()
    blobs  = [b for b in client.list_blobs(GCS_BUCKET, prefix=prefix)
              if b.name.endswith(".tsv")]
    if not blobs:
        raise FileNotFoundError(f"No TSV under gs://{GCS_BUCKET}/{prefix}")
    log.info("Loading: gs://%s/%s", GCS_BUCKET, blobs[0].name)
    return blobs[0].download_as_bytes()


# -- Column helpers ---------------------------------------------------------

def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = (
        df.columns
        .str.lstrip("#")
        .str.strip()
        .str.lower()
        .str.replace(r"[() /%]+", "_", regex=True)
        .str.rstrip("_")
    )
    return df


def extract_country(geo: str) -> str:
    if not geo or pd.isna(geo):
        return "Unknown"
    return str(geo).split(":")[0].strip()


def extract_year(date_str: str) -> int | None:
    if not date_str or pd.isna(date_str):
        return None
    m = re.search(r"\b(19|20)\d{2}\b", str(date_str))
    return int(m.group()) if m else None


# -- Isolate DataFrame ------------------------------------------------------
# Each tuple: (output_column_name, [candidate_input_names_in_priority_order])
# The first candidate found in the actual file wins.

ISOLATE_COLUMN_MAP = [
    ("biosample_acc",    ["biosample_acc", "accession", "sample_name"]),
    ("organism",         ["scientific_name", "organism", "taxgroup_name"]),
    ("strain",           ["strain", "isolate_name", "isolate"]),
    ("collection_date",  ["collection_date", "collection_year"]),
    ("geo_loc_name",     ["geo_loc_name", "geo_loc_name_country_region", "geographic_location"]),
    ("isolation_source", ["isolation_source", "source_type"]),
    ("host",             ["host", "host_disease"]),
    ("serovar",          ["serovar", "serotype", "antigen_formula"]),
    ("amr_genotypes",    ["amr_genotypes", "amr_genotypes_core"]),
]


def build_isolates_df(df: pd.DataFrame) -> pd.DataFrame:
    out = {}
    for out_col, candidates in ISOLATE_COLUMN_MAP:
        for c in candidates:
            if c in df.columns:
                out[out_col] = df[c]
                break
        else:
            out[out_col] = pd.Series([None] * len(df), index=df.index)

    result = pd.DataFrame(out)
    result["country"]         = result["geo_loc_name"].apply(extract_country)
    result["collection_year"] = result["collection_date"].apply(extract_year)

    result = result.dropna(subset=["biosample_acc"])
    result = result[result["biosample_acc"].str.strip() != ""]
    result = result.drop_duplicates(subset=["biosample_acc"])

    log.info("Isolates: %d rows | %d countries | years %s-%s",
             len(result),
             result["country"].nunique(),
             result["collection_year"].min(),
             result["collection_year"].max())
    return result


# -- Gene DataFrame ---------------------------------------------------------

# (regex_prefix, (drug_class, subclass))
# Pseudomonas-relevant genes are listed first; broad patterns last.
GENE_CLASS_MAP = [
    (r"^blaOXA",              ("BETA-LACTAM",    "carbapenemase / OXA")),
    (r"^blaPDC",              ("BETA-LACTAM",    "AmpC cephalosporinase")),
    (r"^blaVEB",              ("BETA-LACTAM",    "ESBL")),
    (r"^mex[A-Z]",            ("EFFLUX",         "MexAB-OprM / RND pump")),
    (r"^opr[A-Z]",            ("EFFLUX",         "outer membrane porin")),
    (r"^fosA",                ("FOSFOMYCIN",     "fosfomycin")),
    (r"^catB",                ("PHENICOL",       "chloramphenicol")),
    (r"^crpP",                ("QUINOLONE",      "ciprofloxacin")),
    (r"^rmtB|^rmtC|^rmtD",   ("AMINOGLYCOSIDE", "16S methyltransferase")),
    (r"^bla",                 ("BETA-LACTAM",    "beta-lactam")),
    (r"^NDM",                 ("BETA-LACTAM",    "carbapenemase / NDM")),
    (r"^KPC",                 ("BETA-LACTAM",    "carbapenemase / KPC")),
    (r"^VIM",                 ("BETA-LACTAM",    "carbapenemase / VIM")),
    (r"^IMP",                 ("BETA-LACTAM",    "carbapenemase / IMP")),
    (r"^CTX-M",               ("BETA-LACTAM",    "ESBL")),
    (r"^tet|^Tet",            ("TETRACYCLINE",   "tetracycline")),
    (r"^aac|^aph|^ant|^aad", ("AMINOGLYCOSIDE", "aminoglycoside")),
    (r"^sul",                 ("SULFONAMIDE",    "sulfonamide")),
    (r"^dfrA",                ("TRIMETHOPRIM",   "trimethoprim")),
    (r"^qnr|^oqx",            ("QUINOLONE",      "quinolone")),
    (r"^mph|^erm",            ("MACROLIDE",      "macrolide")),
    (r"^cat",                 ("PHENICOL",       "chloramphenicol")),
    (r"^mcr",                 ("COLISTIN",       "colistin")),
    (r"^van",                 ("GLYCOPEPTIDE",   "vancomycin")),
]


def classify_gene(gene: str) -> tuple[str, str]:
    for pattern, result in GENE_CLASS_MAP:
        if re.match(pattern, gene, re.IGNORECASE):
            return result
    return ("OTHER", "unclassified")


def build_genes_df(isolates_df: pd.DataFrame) -> pd.DataFrame:
    """Explode amr_genotypes string into one row per gene per isolate."""
    records = []
    for _, row in isolates_df.iterrows():
        geno = row.get("amr_genotypes", "")
        if not geno or pd.isna(geno):
            continue
        genes = [g.strip().split("=")[0] for g in str(geno).split(",") if g.strip()]
        for gene in genes:
            drug_class, subclass = classify_gene(gene)
            records.append({
                "biosample_acc": row["biosample_acc"],
                "gene_symbol":   gene,
                "drug_class":    drug_class,
                "subclass":      subclass,
            })

    genes_df = pd.DataFrame(records)
    if not genes_df.empty:
        log.info("Genes: %d rows | %d unique | %d isolates with hits",
                 len(genes_df),
                 genes_df["gene_symbol"].nunique(),
                 genes_df["biosample_acc"].nunique())
    else:
        log.warning("No gene rows produced — check amr_genotypes column content")
    return genes_df


# -- Main ------------------------------------------------------------------

def main() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw       = load_raw()
    df        = clean_columns(pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str))
    log.info("Raw: %s | cols: %s", df.shape, df.columns.tolist())

    isolates_df = build_isolates_df(df)
    genes_df    = build_genes_df(isolates_df)
    return isolates_df, genes_df


if __name__ == "__main__":
    iso, genes = main()
    print("\n-- Isolates sample --")
    print(iso.head(3).to_string())
    print("\n-- Genes sample --")
    print(genes.head(8).to_string())


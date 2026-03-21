print("Script started")

"""
fetch_amr.py
-----------------

This script downloads the latest AMR metadata file from the
NCBI Pathogen Detection FTP server for a selected organism.

Steps performed:
1. Connect to the NCBI FTP server
2. Locate latest AMR metadata file
3. Download file
4. Save raw copy locally
5. Upload raw copy to GCS (raw/{ORGANISM}/{timestamp}/)
6. Load into pandas
7. Clean columns
8. Basic validation

Environment variables:
    ORGANISM    -- NCBI organism folder (default: Pseudomonas_aeruginosa)
    GCS_BUCKET  -- GCS bucket name for archival
    GOOGLE_APPLICATION_CREDENTIALS -- path to GCP service-account JSON

Run:
    python fetch_amr.py
"""

import os
import io
import re
import ftplib
import logging
from datetime import datetime, UTC

import pandas as pd
import google.auth
from google.cloud import storage


# -----------------------------
# Logging configuration
# -----------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)

log = logging.getLogger(__name__)


# -----------------------------
# Configuration
# -----------------------------
NCBI_FTP_HOST = "ftp.ncbi.nlm.nih.gov"

ORGANISM = os.getenv("ORGANISM", "Pseudomonas_aeruginosa")
GCS_BUCKET = os.getenv("GCS_BUCKET", "amr-data-pipeline-data-lake")

FTP_AMR_PATH = f"/pathogen/Results/{ORGANISM}/latest_snps/AMR"

RAW_DIR = "data"
os.makedirs(RAW_DIR, exist_ok=True)


# -----------------------------
# FTP Connection
# -----------------------------
def connect_ftp():
    log.info("Connecting to NCBI FTP server...")
    ftp = ftplib.FTP(NCBI_FTP_HOST)
    ftp.login()
    log.info("Connected successfully")
    return ftp


# -----------------------------
# Find latest file
# -----------------------------
def find_latest_file(ftp, path, pattern):
    log.info("Navigating to %s", path)
    ftp.cwd(path)

    files = ftp.nlst()

    matched = [f for f in files if re.search(pattern, f)]

    if not matched:
        raise Exception(f"No matching files found in {path}")

    matched.sort(reverse=True)
    latest = matched[0]

    log.info("Latest file found: %s", latest)
    return latest


# -----------------------------
# Download file
# -----------------------------
def download_file(ftp, filename):
    log.info("Downloading %s...", filename)

    buffer = io.BytesIO()
    ftp.retrbinary(f"RETR {filename}", buffer.write)
    buffer.seek(0)

    size_mb = buffer.getbuffer().nbytes / 1_048_576
    log.info("Download complete (%.2f MB)", size_mb)

    return buffer.getvalue()


# -----------------------------
# Clean column names
# -----------------------------
def clean_columns(df):
    df.columns = (
        df.columns
        .str.lstrip("#")
        .str.strip()
        .str.lower()
        .str.replace(r"[() /]+", "_", regex=True)
    )
    return df


# -----------------------------
# GCS Upload
# -----------------------------
def upload_to_gcs(data: bytes, gcs_path: str, run_timestamp: str):
    # Check authentication awareness
    if not os.getenv("GOOGLE_APPLICATION_CREDENTIALS"):
        log.warning(
            "GOOGLE_APPLICATION_CREDENTIALS not set. Ensure ADC is configured."
        )

    credentials, project = google.auth.default()
    client = storage.Client(credentials=credentials, project=project)

    bucket = client.bucket(GCS_BUCKET)
    if not bucket.exists():
        raise ValueError(
            f"GCS bucket '{GCS_BUCKET}' does not exist or is not accessible"
        )

    blob = bucket.blob(gcs_path)

    # Add metadata for traceability
    blob.metadata = {
        "organism": ORGANISM,
        "timestamp": run_timestamp,
        "source": "ncbi_pathogen_detection",
    }

    blob.upload_from_string(
        data,
        content_type="text/tab-separated-values; charset=utf-8",
    )

    log.info("Uploaded to gs://%s/%s", GCS_BUCKET, gcs_path)


# -----------------------------
# Main pipeline
# -----------------------------
def main():
    run_timestamp = datetime.now(UTC).strftime("%Y%m%d_%H%M%S")

    ftp = connect_ftp()

    # -----------------------------
    # AMR FILE
    # -----------------------------
    amr_file = find_latest_file(
        ftp,
        FTP_AMR_PATH,
        r".*amr.*\.tsv",
    )

    amr_bytes = download_file(ftp, amr_file)

    ftp.quit()

    # -----------------------------
    # SAVE RAW FILE LOCALLY
    # -----------------------------
    amr_raw_path = f"{RAW_DIR}/{run_timestamp}_{amr_file}"

    with open(amr_raw_path, "wb") as f:
        f.write(amr_bytes)

    log.info("Saved AMR raw file locally")

    # -----------------------------
    # UPLOAD TO GCS
    # -----------------------------
    gcs_path = f"raw/{ORGANISM}/{run_timestamp}/{amr_file}"

    try:
        upload_to_gcs(amr_bytes, gcs_path, run_timestamp)
    except Exception as e:
        log.error("Failed to upload to GCS: %s", e)
        raise

    # -----------------------------
    # LOAD DATA
    # -----------------------------
    df_amr = pd.read_csv(io.BytesIO(amr_bytes), sep="\t", dtype=str)
    df_amr = clean_columns(df_amr)

    log.info("AMR shape: %s", df_amr.shape)

    # -----------------------------
    # BASIC VALIDATION
    # -----------------------------
    required_cols = ["label", "geo_loc_name", "amr_genotypes"]

    missing = [col for col in required_cols if col not in df_amr.columns]

    if missing:
        log.warning(f"Missing expected columns: {missing}")
    else:
        log.info("All key columns present")

    log.info("Fetch complete. GCS path: gs://%s/%s", GCS_BUCKET, gcs_path)
    log.info("Pipeline completed successfully.")

    return df_amr, gcs_path


if __name__ == "__main__":
    main()
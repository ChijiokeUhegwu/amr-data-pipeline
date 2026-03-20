print("Script started")

"""
ingestion.py
-----------------

This script downloads the latest AMR metadata file from the
NCBI Pathogen Detection FTP server for a selected organism.

Steps performed:
1. Connect to the NCBI FTP server
2. Locate latest AMR metadata file
3. Download file
4. Save raw copy
5. Load into pandas
6. Clean columns

Run:
    python ingestion.py
"""

import os
import io
import re
import ftplib
import logging
from datetime import datetime, UTC

import pandas as pd


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
        r".*amr.*\.tsv"
    )

    amr_bytes = download_file(ftp, amr_file)

    ftp.quit()

    # -----------------------------
    # SAVE RAW FILES
    # -----------------------------
    amr_raw_path = f"{RAW_DIR}/{run_timestamp}_{amr_file}"

    with open(amr_raw_path, "wb") as f:
        f.write(amr_bytes)

    log.info("Saved AMR raw files")

    # -----------------------------
    # LOAD AMR DATA
    # -----------------------------
    df_amr = pd.read_csv(io.BytesIO(amr_bytes), sep="\t", dtype=str)
    df_amr = clean_columns(df_amr)

    log.info("AMR shape: %s", df_amr.shape)
    log.info("Pipeline completed successfully.")


if __name__ == "__main__":
    main()
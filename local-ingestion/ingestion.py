print("Script started")
"""
fetch_ncbi_amr.py
-----------------

This script downloads the latest metadata file from the
NCBI Pathogen Detection FTP server for a selected organism.

Steps performed:
1. Connect to the NCBI FTP server
2. Locate the latest metadata TSV file
3. Download the file
4. Load it into pandas
5. Filter isolates by country
6. Save raw and filtered copies locally

This version is designed for LOCAL testing before deploying
to GCP.

Run:
    python fetch_ncbi_amr.py
"""

import os
import io
import re
import ftplib
import logging
from datetime import datetime

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

# Organism can be changed later
ORGANISM = os.getenv("ORGANISM", "Pseudomonas_aeruginosa")

# FTP location of metadata
FTP_BASE_PATH = f"/pathogen/Results/{ORGANISM}/latest_snps/Metadata"


# Local directories
RAW_DIR = "data/raw"
FILTERED_DIR = "data/filtered"

os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(FILTERED_DIR, exist_ok=True)


# Countries we want to retain
AFRICAN_COUNTRIES = [
    "Nigeria","South Africa","Kenya","Ethiopia","Egypt",
    "Ghana","Tanzania","Uganda","Cameroon","Senegal",
]

EUROPEAN_COUNTRIES = [
    "United Kingdom","Germany","France","Netherlands",
    "Sweden","Denmark","Norway","Spain","Italy"
]

TARGET_COUNTRIES = set(AFRICAN_COUNTRIES + EUROPEAN_COUNTRIES)


# ---------------------------------------------------
# FTP Connection
# ---------------------------------------------------

def connect_ftp():
    """
    Connect to NCBI FTP server using anonymous login.
    """
    log.info("Connecting to NCBI FTP server...")
    ftp = ftplib.FTP(NCBI_FTP_HOST)
    ftp.login()
    log.info("Connected successfully")
    return ftp


# ---------------------------------------------------
# Find metadata file
# ---------------------------------------------------

def find_metadata_file(ftp):
    """
    Navigate to the metadata directory and find the
    latest metadata TSV file.
    """

    log.info("Navigating to %s", FTP_BASE_PATH)
    ftp.cwd(FTP_BASE_PATH)

    files = ftp.nlst()

    # Find metadata files
    metadata_files = [f for f in files if re.search(r"PDG.*metadata\.tsv", f)]

    if not metadata_files:
        raise Exception("No metadata files found.")

    metadata_files.sort(reverse=True)

    latest_file = metadata_files[0]

    log.info("Latest metadata file found: %s", latest_file)

    return latest_file


# ---------------------------------------------------
# Download file
# ---------------------------------------------------

def download_file(ftp, filename):
    """
    Download the file from FTP into memory.
    """

    log.info("Downloading file...")

    buffer = io.BytesIO()

    ftp.retrbinary(f"RETR {filename}", buffer.write)

    buffer.seek(0)

    size_mb = buffer.getbuffer().nbytes / 1_048_576

    log.info("Download complete (%.2f MB)", size_mb)

    return buffer.getvalue()


# ---------------------------------------------------
# Extract country
# ---------------------------------------------------

def extract_country(location):
    """
    Extract the country name from geo_loc_name.

    Examples:
    Nigeria
    Nigeria: Lagos
    Nigeria: Lagos, Ikoyi
    """

    if pd.isna(location):
        return "Unknown"

    return location.split(":")[0].strip()


# ---------------------------------------------------
# Filter isolates by country
# ---------------------------------------------------

def filter_geography(df):
    """
    Filter dataframe to selected countries.
    """

    df["country"] = df["geo_loc_name"].apply(extract_country)

    filtered = df[df["country"].isin(TARGET_COUNTRIES)].copy()

    log.info(
        "Rows before filtering: %d | after filtering: %d",
        len(df),
        len(filtered),
    )

    return filtered


# ---------------------------------------------------
# Main pipeline
# ---------------------------------------------------

def main():

    run_timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")

    # 1. Connect to FTP
    ftp = connect_ftp()

    # 2. Find metadata file
    filename = find_metadata_file(ftp)

    # 3. Download file
    raw_bytes = download_file(ftp, filename)

    ftp.quit()

    # 4. Save raw file locally
    raw_path = f"{RAW_DIR}/{run_timestamp}_{filename}"

    with open(raw_path, "wb") as f:
        f.write(raw_bytes)

    log.info("Raw file saved to %s", raw_path)

    # 5. Load into pandas
    df = pd.read_csv(
        io.BytesIO(raw_bytes),
        sep="\t",
        low_memory=False,
        dtype=str
    )

    # Normalize column names
    df.columns = (
        df.columns
        .str.lstrip("#")
        .str.strip()
        .str.lower()
        .str.replace(r"[() /]+", "_", regex=True)
    )

    log.info("Dataframe shape: %s", df.shape)

    # 6. Filter geography
    filtered_df = filter_geography(df)

    # 7. Save filtered dataset
    filtered_path = f"{FILTERED_DIR}/{run_timestamp}_filtered_metadata.tsv"

    filtered_df.to_csv(filtered_path, sep="\t", index=False)

    log.info("Filtered data saved to %s", filtered_path)

    log.info("Pipeline completed successfully.")


if __name__ == "__main__":
    main()
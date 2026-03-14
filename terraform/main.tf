terraform {
  required_providers {
    google = {
      source  = "hashicorp/google"
      version = "~> 5.0"
    }
  }
}

provider "google" {
  project = var.project_id
  region  = var.region
}

resource "google_storage_bucket" "amr_raw_bucket" {
  name          = "${var.project_id}-data-lake"
  location      = var.region
  force_destroy = true

  uniform_bucket_level_access = true
}

resource "google_bigquery_dataset" "data_warehouse" {
  dataset_id = "data_warehouse"
  location   = var.region
}

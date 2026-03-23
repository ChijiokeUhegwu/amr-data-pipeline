# base Docker image that we will build on
FROM python:3.12-slim

# Copy uv binary from official uv image (multi-stage build pattern)
COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/

# set up the working directory inside the container
WORKDIR /app

# Add virtual environment to PATH so we can use installed packages
ENV PATH="/app/.venv/bin:$PATH"

# Copy dependency files first (better layer caching)
COPY "pyproject.toml" "uv.lock" ".python-version" ./

# Install dependencies from lock file (ensures reproducible builds)
RUN uv sync --locked

# Copy pipeline scripts
COPY ingestion/fetch_amr.py   .
COPY ingestion/transform.py   .
COPY ingestion/load_to_db.py  .

# Default entrypoint runs the full pipeline.
ENTRYPOINT ["python", "/app/load_to_db.py"]

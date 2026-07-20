# Hosting the DepMap Dependency Explorer

This is a Flask app, not a static site. Host it as a Python web service so the
server can keep `OPENAI_API_KEY` private, read the raw DepMap matrices, and
persist custom stratifiers.

## Runtime

Use Gunicorn in production:

```bash
bash scripts/start.sh
```

Minimum environment:

```bash
SECRET_KEY=<long random value>
OPENAI_API_KEY=<OpenAI API key>
OPENAI_MODEL=gpt-5.5
ALLOW_SELF_SIGNUP=false
ADMIN_EMAIL=<initial admin email>
ADMIN_PASSWORD=<initial admin password>
ADMIN_FIRST_NAME=<initial admin first name>
DB_PATH=/var/data/app.db
DEPMAP_DATA_DIR=/var/data/depmap
STRATIFIER_DATASET_DIR=/var/data/stratifier_sources
MAX_STRATIFIER_DATASET_BYTES=26214400
MAX_STRATIFIER_DATASET_ROWS=250000
```

`DB_PATH` must live on persistent storage if you want saved stratifiers and
users to survive redeploys. `DEPMAP_DATA_DIR` must contain the raw DepMap CSVs
used to compute new stratifiers:

- `CRISPRGeneEffect.csv`
- `D2_combined_gene_dep_scores.csv`
- `Model.csv`

`STRATIFIER_DATASET_DIR` stores the source tables selected by ChatGPT web
search. Keep it on persistent storage so saved provenance remains auditable
across redeploys. Downloads are restricted to public HTTP(S) addresses and to
the configured byte and row limits.

## Populating DepMap Data

On a host with `DEPMAP_DATA_DIR` pointing at a persistent disk, run:

```bash
python3 scripts/build_hpv_dependency_data.py
```

The downloader will populate the raw CSV cache and refresh the static summary
used by the built-in dropdown analyses.

## Render / Railway Shape

For Render, Railway, Fly.io, or a VPS, use:

```text
Build: pip install -r requirements.txt
Start: bash scripts/start.sh
```

Then configure the environment variables above. Make sure `DB_PATH` and
`DEPMAP_DATA_DIR` point at a persistent disk or volume, not an ephemeral
checkout directory.

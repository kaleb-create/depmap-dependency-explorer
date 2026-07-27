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
AUTH_ENABLED=false
ALLOW_SELF_SIGNUP=false
ADMIN_EMAIL=<initial admin email>
ADMIN_PASSWORD=<initial admin password>
ADMIN_FIRST_NAME=<initial admin first name>
DB_PATH=/var/data/app.db
DATABASE_URL=<Render Postgres internal database URL>
DEPMAP_DATA_DIR=/var/data/depmap
STRATIFIER_DATASET_DIR=/var/data/stratifier_sources
MAX_STRATIFIER_DATASET_BYTES=26214400
MAX_STRATIFIER_DATASET_ROWS=250000
```

`AUTH_ENABLED=false` runs the explorer without login or registration. Set it to
`true` later to restore session-based access control and user administration.

For hosted deployments, `DATABASE_URL` is the preferred account and stratifier
store. Render can populate it from a linked Postgres database. When
`DATABASE_URL` is set, the app uses Postgres and ignores `DB_PATH`.

Without `DATABASE_URL`, `DB_PATH` must live on persistent storage if you want
saved stratifiers and users to survive redeploys. Render's ordinary web-service
filesystem is ephemeral, so leaving `DB_PATH` at the repository default will
erase registered accounts on a redeploy or restart. `DEPMAP_DATA_DIR` must
contain the raw DepMap CSVs used to compute new stratifiers:

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

Then configure the environment variables above. Connect a Render Postgres
database and expose its internal URL as `DATABASE_URL`. Alternatively, attach a
persistent disk and point `DB_PATH` and `DEPMAP_DATA_DIR` into its mount path,
not an ephemeral checkout directory.

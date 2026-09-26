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
GUNICORN_TIMEOUT=300
DEPMAP_PROVISION_ON_START=true
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

When `DB_PATH` is omitted, an existing `/var/data` disk is now preferred
automatically. Otherwise an explicitly configured `DEPMAP_DATA_DIR` also
provides the default database directory. An existing repository-local database
is migrated only when the destination database does not yet exist. Explicit
`DB_PATH` and `DATABASE_URL` settings always take precedence. Attaching a disk
still requires Render account access; creating a directory is not a substitute
for persistent storage.

- `CRISPRGeneEffect.csv`
- `D2_combined_gene_dep_scores.csv`
- `Model.csv`

On startup, `wsgi.py` checks for the three raw DepMap files and downloads any
missing inputs from immutable official Figshare release URLs before Gunicorn
accepts traffic. The runtime inputs are DepMap Public 24Q4 Model/CRISPR data and
DEMETER2 Data v6 siRNA data, with fixed checksum validation. Without a
persistent disk this cache is rebuilt after every deploy or restart. A disk
mounted at `/var/data` with `DEPMAP_DATA_DIR=/var/data/depmap` avoids repeat
downloads.

`STRATIFIER_DATASET_DIR` stores the source tables selected by ChatGPT web
search. Keep it on persistent storage so saved provenance remains auditable
across redeploys. Downloads are restricted to public HTTP(S) addresses and to
the configured byte and row limits.

## Populating DepMap Data

Custom stratifier requests are saved to `stratifier_jobs` in the configured
database and return immediately. The web process runs one background build
at a time; a file lock in `DEPMAP_DATA_DIR` coordinates Gunicorn workers on
the same host. Reloading the Stratifiers page resumes unfinished jobs after a
restart. Three repeated server interruptions mark a job failed with a retry
option. This design targets a single Render instance with persistent storage.

The planner can use Model.csv, the bundled `stratifier_catalog.json` (including
classified negative cohorts and source definitions), or an external CSV/TSV/JSON
cohort table. Mapping failures are sent back for up to two corrections. Unknown
or conflicting external calls are excluded, and at least three observations
and 50% coverage per cohort are required for every plotted gene. A saved result
includes the actual source and dependency releases; cohort frequencies describe
the classified cell lines, not patient prevalence.

The page displays build progress and retains failed prompts for retry. Completed
jobs link directly to their saved analysis. Custom matrices load on demand from
`/api/dependency-analysis/custom-<id>`.

Runtime provisioning only downloads the three inputs needed for custom
stratifiers. To refresh every built-in analysis from all upstream sources, run:

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

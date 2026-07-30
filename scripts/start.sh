#!/usr/bin/env bash
set -euo pipefail

exec gunicorn --timeout "${GUNICORN_TIMEOUT:-300}" --bind "0.0.0.0:${PORT:-8000}" wsgi:app

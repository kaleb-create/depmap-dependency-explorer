#!/usr/bin/env bash
set -euo pipefail

exec gunicorn --bind "0.0.0.0:${PORT:-8000}" wsgi:app

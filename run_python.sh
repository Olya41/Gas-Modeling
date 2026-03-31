#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p output/data output/plots

PYTHONPATH=python python3 python/main.py
python3 analysis/check.py

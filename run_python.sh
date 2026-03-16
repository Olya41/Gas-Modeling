#!/bin/bash
set -e
cd "$(dirname "$0")"

PYTHONPATH=python python3 python/main.py && python3 analysis/check.py

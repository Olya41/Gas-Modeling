#!/bin/bash
set -e
cd "$(dirname "$0")"

[ -d build ] || cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -Wno-dev > /dev/null
cmake --build build -j$(nproc)

./build/md_simulation && python3 analysis/check.py

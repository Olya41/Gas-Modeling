"""Корень репозитория и чтение params (общее для free_run, msd_plots, check)."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def read_param(name, default=None):
    with open(ROOT / "params") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                if k.strip() == name:
                    return v.strip()
    return default

from pathlib import Path
import subprocess
import sys

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
BLACK_EXCLUDES = r"(^|/)(env|\.venv|venv|build|dist|wheelhouse)(/|$)|(^|/)JupiterMag/__data/libjupitermag(/|$)"
FLAKE8_EXCLUDES = "env,.venv,venv,build,dist,wheelhouse,JupiterMag/__data/libjupitermag"


def _run_check(*args):
    result = subprocess.run(
        [sys.executable, "-m", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        pytest.fail(
            f"Command failed: {sys.executable} -m {' '.join(args)}\n"
            f"Exit code: {result.returncode}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def test_black_check_repo():
    _run_check(
        "black",
        "--line-length",
        "160",
        "--extend-exclude",
        BLACK_EXCLUDES,
        "--check",
        ".",
    )


def test_flake8_repo():
    _run_check("flake8", "--exclude", FLAKE8_EXCLUDES, ".")

#!/usr/bin/env python3
"""
Prepare a local KOspec pipeline checkout for first use.

This script is intentionally small and cross-platform so observers can run it
right after cloning the GitHub repository on an observation PC.
"""

import argparse
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent

DIRECTORIES = [
    ROOT / "data" / "spectra",
    ROOT / "outputs" / "quicklook",
    ROOT / "outputs" / "flux_calibrated",
    ROOT / "logs",
]

REQUIRED_MODULES = {
    "numpy": "numpy",
    "scipy": "scipy",
    "astropy": "astropy",
    "matplotlib": "matplotlib",
}


def module_exists(module_name):
    return importlib.util.find_spec(module_name) is not None


def create_directories():
    for directory in DIRECTORIES:
        directory.mkdir(parents=True, exist_ok=True)
        gitkeep = directory / ".gitkeep"
        if not gitkeep.exists():
            gitkeep.touch()


def create_config():
    config_path = ROOT / "config.yaml"
    example_path = ROOT / "config.example.yaml"
    if config_path.exists() or not example_path.exists():
        return False

    shutil.copyfile(example_path, config_path)
    return True


def missing_modules():
    return [
        package
        for package, module_name in REQUIRED_MODULES.items()
        if not module_exists(module_name)
    ]


def install_dependencies():
    requirements = ROOT / "requirements.txt"
    subprocess.check_call([
        sys.executable,
        "-m",
        "pip",
        "install",
        "-r",
        str(requirements),
    ])


def print_next_steps(config_created, missing):
    print("\nKOspec local setup complete.")
    print(f"\nRepository:")
    print(f"  {ROOT}")
    print("\nInput FITS directory:")
    print("  data/spectra/")
    print("\nOutput directories:")
    print("  outputs/quicklook/")
    print("  outputs/flux_calibrated/")

    if config_created:
        print("\nCreated local config:")
        print("  config.yaml")

    if missing:
        print("\nMissing Python packages:")
        print("  " + ", ".join(missing))
        print("\nInstall them with:")
        print("  python3 -m pip install -r requirements.txt")

    print("\nRun quicklook:")
    print("  python3 main.py")
    print("\nRun quicklook + flux calibration:")
    print("  python3 main_all.py")
    print("\nUse an external observation data directory:")
    print("  python3 main_all.py --spectra-dir /path/to/fits")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare this KOspec checkout for local execution."
    )
    parser.add_argument(
        "--install-deps",
        action="store_true",
        help="Install Python dependencies with pip after creating directories.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if sys.version_info < (3, 8):
        print("Python 3.8 or newer is recommended.", file=sys.stderr)

    create_directories()
    config_created = create_config()

    missing = missing_modules()
    if missing and args.install_deps:
        install_dependencies()
        missing = missing_modules()

    print_next_steps(config_created, missing)
    return 0 if not missing else 1


if __name__ == "__main__":
    sys.exit(main())

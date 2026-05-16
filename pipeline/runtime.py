"""
Runtime environment helpers for command-line execution.
"""

import os
import warnings
from pathlib import Path


def configure_matplotlib_cache():
    """Use a local Matplotlib cache when the user home is not writable."""
    warnings.filterwarnings(
        "ignore",
        message="Unable to import Axes3D.*",
        category=UserWarning,
        module="matplotlib.projections",
    )

    if os.environ.get("MPLCONFIGDIR"):
        return

    root = Path(__file__).resolve().parents[1]
    cache_dir = root / "logs" / "matplotlib"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(cache_dir)

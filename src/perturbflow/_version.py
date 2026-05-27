"""Single source of truth for the package version.

Kept as a plain string (not derived from pyproject at runtime) so that
``perturbflow.__version__`` works in zipapps and Snakemake-rendered Dockers
where ``importlib.metadata`` can be slow or unavailable.
"""

from __future__ import annotations

__version__ = "0.1.0"

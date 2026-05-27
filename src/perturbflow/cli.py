"""``perturbflow`` command-line interface.

Subcommands
-----------
- ``perturbflow run --config CONFIG`` — end-to-end pipeline.
- ``perturbflow validate-config --config CONFIG`` — parse without running.
- ``perturbflow version`` — print version and exit.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

from perturbflow._version import __version__


def _configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stderr,
    )


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="perturbflow")
def main() -> None:
    """PerturbFlow: opinionated Perturb-seq analysis pipeline."""


@main.command("run")
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the run YAML config.",
)
@click.option(
    "--outdir",
    default=None,
    type=click.Path(path_type=Path),
    help="Override config.run.outdir.",
)
@click.option(
    "--log-level",
    default=None,
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    help="Override log level from config.",
)
def run_cmd(config_path: Path, outdir: Path | None, log_level: str | None) -> None:
    """Run the full PerturbFlow pipeline."""
    from perturbflow.config import load_config
    from perturbflow.pipeline import run as run_pipeline

    cfg = load_config(config_path)
    _configure_logging(log_level or cfg.run.log_level)
    artifacts = run_pipeline(cfg, config_path=config_path, outdir=outdir)
    out = outdir or cfg.run.outdir
    click.echo(f"PerturbFlow complete. Outputs in: {out}")
    if artifacts.report_path:
        click.echo(f"Report: {artifacts.report_path}")


@main.command("validate-config")
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
def validate_config_cmd(config_path: Path) -> None:
    """Parse a config file and exit (non-zero on error)."""
    from perturbflow.config import load_config

    load_config(config_path)
    click.echo(f"OK: {config_path}")


@main.command("version")
def version_cmd() -> None:
    """Print the package version."""
    click.echo(__version__)


if __name__ == "__main__":
    main()

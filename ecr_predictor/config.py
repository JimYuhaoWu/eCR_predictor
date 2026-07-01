"""
Load and validate config.yaml for the refinement pipeline.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

_DEFAULT_CONFIG_PATH = Path(__file__).parents[1] / "config.yaml"


def load_config(path: Path | str | None = None) -> dict[str, Any]:
    """
    Load config.yaml and return as a nested dict.

    Falls back to config.yaml in the repo root if path is None.
    Raises FileNotFoundError if the file does not exist.
    """
    try:
        import yaml
    except ImportError:
        raise ImportError(
            "PyYAML is required for config loading: pip install pyyaml"
        )

    config_path = Path(path) if path else _DEFAULT_CONFIG_PATH
    if not config_path.exists():
        raise FileNotFoundError(
            f"Config file not found: {config_path}\n"
            f"Copy config.yaml from the repo root and edit it."
        )

    with config_path.open(encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    return cfg or {}


def get_af3_config(cfg: dict[str, Any]) -> dict[str, Any]:
    """Return the af3 section of the config, with safe defaults."""
    af3 = cfg.get("af3", {})
    af3.setdefault("backend", "hpcc")
    return af3


def get_fusion_config(cfg: dict[str, Any]) -> dict[str, Any]:
    """Return the fusion section of the config, with safe defaults."""
    fusion = cfg.get("fusion", {})
    fusion.setdefault("orientation", "dbd_n")
    fusion.setdefault("linkers", [])
    fusion.setdefault("self_proteome", "")
    fusion.setdefault("proteome_dir", "data/proteomes")
    fusion.setdefault("junction_window", 5)
    fusion.setdefault("tools", {})
    return fusion

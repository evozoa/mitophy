"""Load and merge YAML configuration; profile overlay; typed accessors."""
from __future__ import annotations
import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import yaml

ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = ROOT / "config"
DATA_DIR = ROOT / "data"


def _deep_merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_yaml(path: Path) -> Any:
    with open(path) as fh:
        return yaml.safe_load(fh)


@dataclass
class Config:
    profile: str = "default"
    root: Path = ROOT
    pipeline: dict = field(default_factory=dict)
    markers: dict = field(default_factory=dict)
    taxa_origin: dict = field(default_factory=dict)
    sampling: dict = field(default_factory=dict)
    supergroups: dict = field(default_factory=dict)
    site: dict = field(default_factory=dict)
    workdir: Path = ROOT / "work"
    resultsdir: Path = ROOT / "results"
    force: bool = False
    threads: int | None = None

    # ---- convenience -------------------------------------------------------
    @property
    def params(self) -> dict:
        """Profile-resolved parameter block (threads, origin, diversification)."""
        prof = self.pipeline.get("profiles", {})
        p = _deep_merge(prof.get("default", {}), prof.get(self.profile, {}) if self.profile != "default" else {})
        if self.threads:
            p["threads"] = self.threads
        return p

    @property
    def tools(self) -> dict:
        return self.pipeline.get("tools", {})

    @property
    def marker_list(self) -> list[dict]:
        return self.markers["markers"]

    def markers_for(self, analysis: str) -> list[dict]:
        return [m for m in self.marker_list if analysis in m.get("use", [])]

    def sub(self, *parts: str) -> Path:
        p = self.workdir.joinpath(*parts)
        p.mkdir(parents=True, exist_ok=True)
        return p


def load_config(profile: str = "default", workdir: str | Path | None = None, resultsdir: str | Path | None = None,
                force: bool = False, threads: int | None = None, config_dir: Path = CONFIG_DIR) -> Config:
    pipeline = load_yaml(config_dir / "pipeline.yml")
    if profile not in pipeline.get("profiles", {}):
        raise SystemExit(f"unknown profile '{profile}'; available: {', '.join(pipeline['profiles'])}")
    wd = Path(workdir) if workdir else ROOT / pipeline.get("workdir", "work")
    rd = Path(resultsdir) if resultsdir else ROOT / pipeline.get("resultsdir", "results")
    return Config(
        profile=profile,
        pipeline=pipeline,
        markers=load_yaml(config_dir / "markers.yml"),
        taxa_origin=load_yaml(config_dir / "taxa_origin.yml"),
        sampling=load_yaml(config_dir / "sampling_diversification.yml"),
        supergroups=load_yaml(DATA_DIR / "taxonomy" / "supergroups.yml"),
        site=load_yaml(config_dir / "site.yml"),
        workdir=wd, resultsdir=rd, force=force, threads=threads,
    )

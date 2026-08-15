"""Subprocess wrappers for external tools with version capture and binary auto-detection."""
from __future__ import annotations
import logging, re, shutil, subprocess, sys
from pathlib import Path

log = logging.getLogger("mitophy")

def which_first(*names: str) -> str | None:
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    return None

def iqtree_bin(cfg_bin: str = "auto") -> str:
    if cfg_bin and cfg_bin != "auto":
        return cfg_bin
    b = which_first("iqtree3", "iqtree2", "iqtree")
    if not b:
        raise SystemExit("IQ-TREE not found (tried iqtree3, iqtree2, iqtree)")
    return b

def fasttree_bin(cfg_bin: str = "auto") -> str:
    if cfg_bin and cfg_bin != "auto":
        return cfg_bin
    b = which_first("FastTreeMP", "fasttree", "FastTree")
    if not b:
        raise SystemExit("FastTree not found")
    return b

def run(cmd: list[str], cwd: Path | None = None, log_file: Path | None = None, check: bool = True, quiet: bool = False,
        stdout_to: Path | None = None) -> subprocess.CompletedProcess:
    cmd = [str(c) for c in cmd]
    log.info("$ %s", " ".join(cmd))
    if stdout_to is not None:
        with open(stdout_to, "w") as out, open(log_file or (str(stdout_to) + ".log"), "w") as err:
            p = subprocess.run(cmd, cwd=cwd, stdout=out, stderr=err, text=True)
    else:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
        if log_file:
            Path(log_file).write_text((p.stdout or "") + "\n" + (p.stderr or ""))
    if check and p.returncode != 0:
        tail = ""
        if stdout_to is None:
            tail = (p.stderr or p.stdout or "")[-2000:]
        elif log_file and Path(log_file).exists():
            tail = Path(log_file).read_text()[-2000:]
        raise RuntimeError(f"command failed ({p.returncode}): {' '.join(cmd)}\n{tail}")
    return p

def _version(cmd: list[str], pattern: str, stderr_ok: bool = True) -> str:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        txt = (p.stdout or "") + ("\n" + p.stderr if stderr_ok else "")
        m = re.search(pattern, txt)
        return m.group(1) if m else txt.strip().splitlines()[0][:80]
    except Exception as e:  # pragma: no cover
        return f"unknown ({e.__class__.__name__})"

def tool_versions(cfg_tools: dict) -> dict:
    v = {"python": sys.version.split()[0]}
    try:
        import Bio, pyhmmer
        v["biopython"] = Bio.__version__; v["pyhmmer"] = pyhmmer.__version__
    except Exception:
        pass
    if which_first("mafft"):
        v["mafft"] = _version(["mafft", "--version"], r"v(\d+\.\d+)")
    if which_first("trimal"):
        v["trimal"] = _version(["trimal", "--version"], r"trimAl v?([\d.a-zA-Z]+)")
    try:
        ib = iqtree_bin(cfg_tools.get("iqtree", {}).get("bin", "auto"))
        v["iqtree"] = Path(ib).name + " " + _version([ib, "--version"], r"IQ-TREE.*?version (\S+)")
    except SystemExit:
        pass
    try:
        fb = fasttree_bin(cfg_tools.get("fasttree", {}).get("bin", "auto"))
        v["fasttree"] = _version([fb], r"FastTree [Vv]ersion (\S+)")
    except SystemExit:
        pass
    if which_first("datasets"):
        v["datasets"] = _version(["datasets", "--version"], r"version:?\s*(\S+)")
    return v

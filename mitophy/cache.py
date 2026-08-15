"""Content hashing and per-stage manifests for idempotent reruns."""
from __future__ import annotations
import hashlib, json, time
from pathlib import Path
from typing import Iterable

def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()

def hash_inputs(files: Iterable[Path] = (), strings: Iterable[str] = (), objs: Iterable[object] = ()) -> str:
    """Combined hash of file contents, strings and JSON-serialisable objects (order-sensitive)."""
    h = hashlib.sha256()
    for f in files:
        f = Path(f)
        h.update(f.name.encode()); h.update(sha256_file(f).encode() if f.exists() else b"MISSING")
    for s in strings:
        h.update(s.encode())
    for o in objs:
        h.update(json.dumps(o, sort_keys=True, default=str).encode())
    return h.hexdigest()

def manifest_path(outdir: Path) -> Path:
    return Path(outdir) / "manifest.json"

def read_manifest(outdir: Path) -> dict:
    p = manifest_path(outdir)
    return json.loads(p.read_text()) if p.exists() else {}

def write_manifest(outdir: Path, input_hash: str, **extra) -> dict:
    m = {"input_hash": input_hash, "finished": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **extra}
    Path(outdir).mkdir(parents=True, exist_ok=True)
    manifest_path(outdir).write_text(json.dumps(m, indent=2, default=str))
    return m

def up_to_date(outdir: Path, input_hash: str, outputs: Iterable[Path] = (), force: bool = False) -> bool:
    if force:
        return False
    m = read_manifest(outdir)
    if m.get("input_hash") != input_hash:
        return False
    return all(Path(o).exists() for o in outputs)

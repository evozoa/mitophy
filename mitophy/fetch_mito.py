"""Stage 1: download RefSeq mitochondrion release files (GenBank flat files), cached by HTTP Last-Modified."""
from __future__ import annotations
import fnmatch, json, logging, re, time
from pathlib import Path
import requests
from .config import Config

log = logging.getLogger("mitophy")

def list_remote(dir_url: str, pattern: str) -> list[str]:
    r = requests.get(dir_url, timeout=60)
    r.raise_for_status()
    names = set(re.findall(r'href="([^"]+)"', r.text))
    return sorted(n for n in names if fnmatch.fnmatch(n, pattern))

def _head(url: str) -> dict:
    r = requests.head(url, timeout=60, allow_redirects=True)
    r.raise_for_status()
    return {"last_modified": r.headers.get("Last-Modified"), "size": int(r.headers.get("Content-Length", 0))}

def _download(url: str, dest: Path, retries: int = 3) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            tmp.replace(dest)
            return
        except Exception as e:  # noqa
            log.warning("download failed (%s/%s): %s", attempt + 1, retries, e)
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"could not download {url}")

def run(cfg: Config) -> dict:
    src = cfg.pipeline["sources"]
    outdir = cfg.sub("mito")
    source_json = outdir / "source.json"
    prev = json.loads(source_json.read_text()) if source_json.exists() else {"files": {}}
    try:
        names = list_remote(src["refseq_mito_dir"], src["refseq_mito_glob"])
    except Exception as e:
        if prev["files"] and all((outdir / n).exists() for n in prev["files"]):
            log.warning("remote listing failed (%s); using cached files", e)
            return prev
        raise
    if not names:
        raise RuntimeError("no RefSeq mitochondrion files matched")
    files = {}
    for n in names:
        url = src["refseq_mito_dir"] + n
        meta = _head(url)
        dest = outdir / n
        cached = prev["files"].get(n, {})
        if dest.exists() and cached.get("last_modified") == meta["last_modified"] and dest.stat().st_size == meta["size"] and not cfg.force:
            log.info("cached: %s (%s)", n, meta["last_modified"])
        else:
            log.info("downloading %s (%.1f MB)", n, meta["size"] / 1e6)
            _download(url, dest)
        files[n] = {**meta, "path": str(dest)}
    # remove stale files from previous releases
    for old in prev["files"]:
        if old not in files and (outdir / old).exists():
            (outdir / old).unlink()
    info = {"files": files, "dir": src["refseq_mito_dir"], "fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    source_json.write_text(json.dumps(info, indent=2))
    return info

def gbff_files(cfg: Config) -> list[Path]:
    source_json = cfg.workdir / "mito" / "source.json"
    if not source_json.exists():
        raise SystemExit("run `mitophy fetch-mito` first")
    info = json.loads(source_json.read_text())
    return [Path(v["path"]) for v in info["files"].values()]

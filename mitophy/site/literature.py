"""Load curated literature trees (Newick + JSON sidecar) and validate them; format references."""
from __future__ import annotations
import html, json, logging
from pathlib import Path
import dendropy

log = logging.getLogger("mitophy")

def load_literature(lit_dir: Path, out_dir: Path) -> dict[str, dict]:
    """Return {tree_id: merged json} and write merged viewer JSON files into out_dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    trees = {}
    for sub in ("origin", "diversification"):
        for nwk in sorted((lit_dir / sub).glob("*.nwk")):
            side = nwk.with_suffix(".json")
            meta = json.loads(side.read_text()) if side.exists() else {}
            t = dendropy.Tree.get(path=str(nwk), schema="newick", suppress_internal_node_taxa=True, preserve_underscores=True)
            tips = {n.taxon.label for n in t.leaf_node_iter()}
            nodes = {n.label for n in t.internal_nodes() if n.label}
            missing = tips - set(meta.get("tips", {}))
            if missing:
                log.warning("%s: tips without metadata: %s", nwk.name, sorted(missing))
            unknown = set(meta.get("nodes", {})) - nodes
            if unknown:
                log.warning("%s: node annotations without matching labels: %s", nwk.name, sorted(unknown))
            data = {**meta, "newick": nwk.read_text().strip(), "n_taxa": len(tips), "kind": "literature", "section": sub}
            (out_dir / f"{nwk.stem}.json").write_text(json.dumps(data, ensure_ascii=False))
            trees[nwk.stem] = data
    return trees

class ReferenceFormatter:
    def __init__(self, refs: dict):
        self.refs = refs

    def short(self, key: str) -> str:
        r = self.refs.get(key)
        if not r:
            return key
        auth = [a.strip() for a in r["authors"].split(",")]
        surname = lambda a: a.split()[0]
        first = surname(auth[0])
        n = len(auth)
        return f"{first} et al. {r['year']}" if n > 2 else (f"{first} & {surname(auth[1])} {r['year']}" if n == 2 else f"{first} {r['year']}")

    def cite(self, key: str) -> str:
        r = self.refs.get(key)
        if not r:
            return f"<span class='badge'>{html.escape(key)}</span>"
        return f"<a href='https://doi.org/{html.escape(r['doi'])}' title='{html.escape(r['title'])}' target='_blank' rel='noopener'>{html.escape(self.short(key))}</a>"

    def html(self, key: str) -> str:
        r = self.refs.get(key)
        if not r:
            return html.escape(key)
        return (f"{html.escape(r['authors'])} ({r['year']}). {html.escape(r['title'])}. <i>{html.escape(r['journal'])}</i> "
                f"{html.escape(str(r.get('volume', '')))}{(': ' + html.escape(str(r['pages']))) if r.get('pages') else ''}. "
                f"<a href='https://doi.org/{html.escape(r['doi'])}' target='_blank' rel='noopener'>doi:{html.escape(r['doi'])}</a>")

"""Lineage parsing and mapping to display groups (data/taxonomy/supergroups.yml)."""
from __future__ import annotations

def lineage_list(lineage: str | list[str]) -> list[str]:
    if isinstance(lineage, list):
        return [x.strip() for x in lineage if x.strip()]
    return [x.strip() for x in lineage.replace("\n", " ").split(";") if x.strip()]

class GroupMapper:
    def __init__(self, supergroups_cfg: dict):
        self.groups = supergroups_cfg["groups"]
        self.by_id = {g["id"]: g for g in self.groups}

    def group_of(self, lineage: str | list[str], organism: str = "") -> dict:
        lin = lineage_list(lineage)
        hay = set(lin) | {organism, organism.split()[0] if organism else ""}
        for g in self.groups:
            for m in g["match"]:
                if m in hay or any(m in x for x in lin):
                    return g
        return {"id": "Unassigned", "color": "#999999", "supergroup": "Unassigned"}

    def color(self, group_id: str) -> str:
        return self.by_id.get(group_id, {}).get("color", "#999999")

    def legend(self, ids: list[str]) -> list[dict]:
        seen = []
        for g in self.groups:
            if g["id"] in ids and g["id"] not in [s["id"] for s in seen]:
                seen.append({"id": g["id"], "color": g["color"], "supergroup": g.get("supergroup", "")})
        return seen

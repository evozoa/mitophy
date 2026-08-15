#!/usr/bin/env python
"""Resolve every DOI in data/literature/references.yml via doi.org (HEAD) and report failures."""
import sys, requests, yaml
refs = yaml.safe_load(open("data/literature/references.yml"))
bad = 0
for k, r in refs.items():
    try:
        resp = requests.head(f"https://doi.org/{r['doi']}", allow_redirects=False, timeout=30, headers={"User-Agent": "mitophy-doi-check"})
        ok = resp.status_code in (301, 302, 303, 307, 308)
    except Exception as e:
        ok = False; resp = None
    print(f"{'ok ' if ok else 'BAD'} {k:24s} {r['doi']}  {getattr(resp,'status_code','-')} {resp.headers.get('Location','')[:80] if resp is not None else ''}")
    bad += not ok
sys.exit(1 if bad else 0)

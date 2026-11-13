# -*- coding: utf-8 -*-
"""33 参考文献 DOI 核验：用 Crossref API 逐条解析 DOI，比对返回标题与本文引用标题。只读。"""
from pathlib import Path
import html
import re
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "reference_check.csv"


def parse_refs():
    t = (ROOT / "reports" / "part4.html").read_text(encoding="utf-8")
    refs = []
    for m in re.finditer(r"<p>\[(\d+)\](.*?)</p>", t, re.S):
        body = html.unescape(re.sub(r"<[^>]+>", "", m.group(2)))
        body = re.sub(r"\s+", " ", body).strip()
        doi = re.search(r"doi:\s*(10\.[^\s<]+?)\.?$", body)
        refs.append((int(m.group(1)), body, doi.group(1).rstrip(".") if doi else ""))
    return refs


def norm(s):
    return {w for w in re.findall(r"[a-z0-9]+", s.lower()) if len(w) > 3}


def check(doi):
    for attempt in range(2):
        try:
            r = requests.get(f"https://api.crossref.org/works/{doi}",
                             timeout=25, headers={"User-Agent": "aidd-audit/1.0 (mailto:user@example.org)"})
            if r.status_code == 404:
                return "DOI_NOT_FOUND", "", "", ""
            if r.status_code != 200:
                time.sleep(1)
                continue
            msg = r.json()["message"]
            title = (msg.get("title") or [""])[0]
            year = ""
            for k in ("published-print", "published-online", "issued"):
                if msg.get(k, {}).get("date-parts"):
                    year = str(msg[k]["date-parts"][0][0])
                    break
            journal = (msg.get("container-title") or [""])[0]
            return "OK", title, year, journal
        except Exception:
            time.sleep(1)
    return "NETWORK_FAIL", "", "", ""


def main():
    rows = []
    for n, body, doi in parse_refs():
        cited_title = body.split(". ")[0] if not doi else body
        if doi:
            status, title, year, journal = check(doi)
            overlap = ""
            if title:
                a, b = norm(title), norm(body)
                overlap = f"{len(a & b)}/{len(a)}"
                if len(a & b) / max(len(a), 1) < 0.5:
                    status = "TITLE_MISMATCH"
        else:
            status, title, year, journal, overlap = "NO_DOI", "", "", "", ""
        rows.append({"ref": n, "doi": doi, "status": status, "crossref_title": title[:90],
                     "crossref_year": year, "crossref_journal": journal[:40],
                     "title_token_overlap": overlap})
        print(f"[{n}] {status:14s} {doi:38s} {title[:60]}")

    import pandas as pd
    pd.DataFrame(rows).to_csv(OUT, index=False, encoding="utf-8-sig")
    bad = [r for r in rows if r["status"] != "OK"]
    print("\n非 OK 条目:", [(r["ref"], r["status"]) for r in bad] or "无")
    print("written", OUT)


if __name__ == "__main__":
    main()
# -*- coding: utf-8 -*-
"""31 附录 A 与参考文献一致性核查。

- 从 logs/model_log.md 解析五折签名基因，计算每个基因的跨折出现次数；
- 与 reports/part4.html 附录 A 表格逐行比对；
- 抽取全文引用编号，检查是否都存在、是否重复、是否缺号。

输出：reports/appendixA_check.csv、reports/citation_check.csv
运行：C:/Anaconda/Anaconda3/python.exe src/31_verify_appendix_refs.py
"""
from pathlib import Path
import re
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LOGS = ROOT / "logs"
REP = ROOT / "reports"


def fold_genes():
    txt = (LOGS / "model_log.md").read_text(encoding="utf-8")
    folds = {}
    for m in re.finditer(r"- fold (\d+)：折内 DEG \d+ 个，签名 \d+ 个.*?基因：(.+)", txt):
        folds[int(m.group(1))] = [g.strip() for g in m.group(2).split(",")]
    return folds


def main():
    folds = fold_genes()
    cnt = {}
    for k, genes in folds.items():
        for g in genes:
            cnt[g] = cnt.get(g, 0) + 1
    df = pd.DataFrame(sorted(cnt.items(), key=lambda x: (-x[1], x[0])),
                      columns=["gene", "n_folds"])
    df["folds"] = df["gene"].map(lambda g: ";".join(str(k) for k, v in sorted(folds.items()) if g in v))
    df.to_csv(REP / "appendixA_check.csv", index=False, encoding="utf-8-sig")
    print("=== 由 model_log 重算的跨折频率 ===")
    print(df.to_string(index=False))

    # 与附录 A 表格比对
    html = (REP / "part4.html").read_text(encoding="utf-8")
    m = re.search(r"<h1[^>]*>Appendix A.*?</table>", html, re.S)
    seg = m.group(0) if m else ""
    rows = re.findall(r"<tr>(.*?)</tr>", seg, re.S)
    appendix, ok = {}, True
    for r in rows:
        cells = [re.sub(r"<[^>]+>", "", c).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", r, re.S)]
        if len(cells) >= 3 and cells[0].isdigit():
            appendix[int(cells[0])] = [g.strip() for g in cells[2].split(",") if g.strip()]
    print("\n=== 附录 A 与日志比对 ===")
    for k in sorted(folds):
        a = appendix.get(k)
        same = a == folds[k]
        ok = ok and same
        print(f"fold {k}: 附录A={a}\n         日志={folds[k]}\n         一致={same}")
    print("ALL_MATCH", ok)

    # 参考文献
    refs = (REP / "part4.html").read_text(encoding="utf-8")
    ref_ids = re.findall(r"<p>\[(\d+)\]", refs)
    used = set()
    for p in ["part1.html", "part2.html", "part3.html", "part4.html"]:
        t = (REP / p).read_text(encoding="utf-8")
        used |= {int(x) for x in re.findall(r"\[(\d+)\]", t)}
    listed = [int(x) for x in ref_ids]
    rows = []
    for i in sorted(set(listed) | used):
        rows.append({"ref": i, "in_list": i in listed, "cited_in_text": i in used,
                     "n_list_occurrences": listed.count(i)})
    cite = pd.DataFrame(rows)
    cite.to_csv(REP / "citation_check.csv", index=False, encoding="utf-8-sig")
    print("\n=== 引用检查：缺号/重复 ===")
    print(cite[(~cite.in_list) | (~cite.cited_in_text) | (cite.n_list_occurrences > 1)].to_string(index=False))
    print(f"listed={len(listed)} unique_listed={len(set(listed))} max={max(listed) if listed else 0}")


if __name__ == "__main__":
    main()
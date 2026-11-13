# -*- coding: utf-8 -*-
"""预处理：log2 核验、探针->基因映射（一基因多探针取行均值最高者）、输出建模就绪表。

输入：data/interim/{GSE}_expr.parquet、data/interim/{GSE}_meta.csv、data/raw/GPL570.annot.gz
输出：data/processed/{GSE}_gene.csv、data/processed/{GSE}_meta.csv、data/processed/{GSE}_probe_map.csv
       logs/preprocess_log.md、logs/notes.md
运行：C:/Anaconda/Anaconda3/python.exe src/03_preprocess.py
说明：两数据集独立处理，不合并、不做 ComBat、不筛基因（除无注释探针）、不做 DEG、不建模。
"""
from pathlib import Path
import gzip
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"
LOGS = ROOT / "logs"
PROCESSED.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)

DATASETS = ["GSE66360", "GSE48060"]
ANNOT = RAW / "GPL570.annot.gz"


def load_gene_symbols(annot_path: Path):
    """从 GPL570.annot 读取 probe_id -> 首个非空 Gene symbol。"""
    with gzip.open(annot_path, "rt", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    start = None
    for i, l in enumerate(lines):
        if l.strip() == "!platform_table_begin":
            start = i + 1
            break
    if start is None:
        raise RuntimeError("未找到 !platform_table_begin")
    header = lines[start].split("\t")
    gi = header.index("Gene symbol") if "Gene symbol" in header else 2
    mapping = {}
    for line in lines[start + 1:]:
        if line.strip() == "!platform_table_end":
            break
        parts = line.split("\t")
        if len(parts) <= gi:
            continue
        pid = parts[0].strip()
        sym = parts[gi].strip()
        first = ""
        for s in sym.split("///"):
            s = s.strip()
            if s:
                first = s
                break
        mapping[pid] = first
    return mapping


def update_notes(notes_path: Path, items):
    existing = notes_path.read_text(encoding="utf-8") if notes_path.exists() else "# 项目记录\n\n"
    if "CEC" in existing and "探针合并规则" in existing:
        return False
    add = "\n## P2 记录\n\n" + "\n".join(f"- {t}" for t in items) + "\n"
    notes_path.write_text(existing + add, encoding="utf-8")
    return True


def main():
    gene_map = load_gene_symbols(ANNOT)
    print(f"[annot] 解析 {len(gene_map)} 个探针注释")
    log = ["# P2 预处理日志", ""]

    for gse in DATASETS:
        print(f"\n===== {gse} =====")
        expr = pd.read_parquet(INTERIM / f"{gse}_expr.parquet")
        meta = pd.read_csv(INTERIM / f"{gse}_meta.csv", encoding="utf-8-sig")
        n_probe_raw, n_sample = expr.shape
        print(f"原始：{n_probe_raw} 探针 x {n_sample} 样本")

        # --- log2 核验（逐样本中位数/99分位/最大值）
        per_median = expr.median(axis=0)
        per_p99 = expr.quantile(0.99, axis=0)
        per_max = expr.max(axis=0)
        overall_median = float(np.nanmedian(expr.values))
        print("逐样本中位数 / 99分位 / 最大值（前 5 个样本示例）：")
        for sid in list(expr.columns)[:5]:
            print(f"  {sid}: median={per_median[sid]:.3f}  p99={per_p99[sid]:.3f}  max={per_max[sid]:.3f}")
        need_log2 = overall_median > 100
        print(f"数据集整体中位数={overall_median:.3f} -> {'未 log2，执行 log2(x+1)' if need_log2 else '已 log2，不变换'}")
        if need_log2:
            expr = np.log2(expr + 1)

        # --- 探针->基因
        gene_col = [gene_map.get(pid, "") for pid in expr.index]
        full = pd.DataFrame({"probe_id": expr.index, "gene_symbol": gene_col})
        full["row_mean"] = expr.mean(axis=1).values
        dropped = int((full["gene_symbol"] == "").sum())
        valid = full[full["gene_symbol"] != ""].copy()
        valid = valid.sort_values("row_mean", ascending=False)
        valid["selected"] = ~valid.duplicated("gene_symbol", keep="first")
        selected_ids = valid[valid["selected"]]["probe_id"].tolist()
        gene_expr = expr.loc[selected_ids].copy()
        gene_expr.index = valid[valid["selected"]]["gene_symbol"].values
        gene_expr.index.name = "gene_symbol"

        full["selected"] = full["probe_id"].isin(selected_ids).astype(int)
        full.to_csv(PROCESSED / f"{gse}_probe_map.csv", index=False, encoding="utf-8-sig")
        gene_expr.to_csv(PROCESSED / f"{gse}_gene.csv", encoding="utf-8-sig")
        meta.to_csv(PROCESSED / f"{gse}_meta.csv", index=False, encoding="utf-8-sig")

        n_gene = gene_expr.shape[0]
        n_merged = n_probe_raw - dropped - n_gene
        loss_rate = (n_probe_raw - n_gene) / n_probe_raw * 100
        print(f"空注释丢弃探针：{dropped}")
        print(f"映射后基因数：{n_gene}（一基因多探针被合并：{n_merged} 个探针）")
        print(f"最终维度：{n_gene} 基因 x {gene_expr.shape[1]} 样本")
        print(f"探针丢失率：{loss_rate:.2f}%")
        print(f"各 role 样本数：{meta['role'].value_counts().to_dict()}")
        print(f"各 group 样本数：{meta['group'].value_counts().to_dict()}")

        log.append(f"## {gse}")
        log.append(f"- 原始探针数：{n_probe_raw}")
        log.append(f"- 空注释丢弃探针：{dropped}")
        log.append(f"- 映射基因数：{n_gene}")
        log.append(f"- 一基因多探针被合并的探针数：{n_merged}")
        log.append(f"- 最终维度：{n_gene} 基因 x {gene_expr.shape[1]} 样本")
        log.append(f"- 探针丢失率：{loss_rate:.2f}%")
        if need_log2:
            log.append(f"- 是否触发二次 log2：是（整体中位数 {overall_median:.3f} > 100，执行 log2(x+1)）")
        else:
            log.append(f"- 是否触发二次 log2：否（整体中位数 {overall_median:.3f} <= 100）")
        log.append(f"- 各 role 样本数：{meta['role'].value_counts().to_dict()}")
        log.append(f"- 各 group 样本数：{meta['group'].value_counts().to_dict()}")
        log.append("")

    log.append("## 探针合并规则变更")
    log.append("- P1 曾记'取 median'；P2 起统一为'取行均值最高的探针'，理由：探针间信号差异主要反映杂交效率，取最强探针可保留该基因的可测信号，与 genefilter::findLargest 一致。")
    log.append("")
    (LOGS / "preprocess_log.md").write_text("\n".join(log), encoding="utf-8")

    notes_items = [
        "GSE66360 样本为 CD146+ 循环内皮细胞（CEC），报告须写 CEC，不得表述为'外周血全血'。",
        "探针合并规则若改为 median，DEG 数量可能变化，未做敏感性分析，写入 Limitations。",
    ]
    update_notes(LOGS / "notes.md", notes_items)

    print(f"\n[03_preprocess] 已写入 logs/preprocess_log.md, notes.md")


if __name__ == "__main__":
    main()
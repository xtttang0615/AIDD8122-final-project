# -*- coding: utf-8 -*-
"""解析 GEO series matrix：表头字段、表达矩阵、样本分组、尺度探测。

输入：data/raw/GSE66360_series_matrix.txt.gz、data/raw/GSE48060_series_matrix.txt.gz
输出：data/interim/<GSE>_expr.parquet、data/interim/<GSE>_meta.csv、logs/parse_log.md
运行：C:/Anaconda/Anaconda3/python.exe src/02_parse_geo.py
说明：本阶段不做 log2、不做探针映射、不做过滤，仅解析与记录。
"""
from pathlib import Path
import gzip
import sys
from collections import Counter, defaultdict

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
LOGS = ROOT / "logs"
INTERIM.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "GSE66360": "GSE66360_series_matrix.txt.gz",
    "GSE48060": "GSE48060_series_matrix.txt.gz",
}

CASE_KEYWORDS = ["acute myocardial infarction", "myocardial infarction", "stemi", "ami", "patient"]
CONTROL_KEYWORDS = ["healthy", "normal", "control"]


def read_lines(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        return f.read().splitlines()


def parse_series_matrix(lines):
    """解析 series matrix，返回 (header_fields: dict[str, list[str]], table_lines: list[str])。"""
    header = {}
    table = []
    in_table = False
    cur = None
    buf = []

    def flush():
        nonlocal buf, cur
        if cur is not None and buf:
            header[cur] = buf
        buf = []
        cur = None

    for line in lines:
        if line.startswith("!"):
            flush()
            if "\t" in line:
                parts = line.split("\t", 1)
                field = parts[0].strip().strip('"')
                inline = parts[1].rstrip("\n")
            else:
                field = line.strip().strip('"')
                inline = None

            if field == "!series_matrix_table_begin":
                in_table = True
                continue
            if field == "!series_matrix_table_end":
                in_table = False
                continue

            cur = field
            if inline is not None:
                buf = [inline]
        else:
            if in_table:
                table.append(line)
            elif cur is not None:
                buf.append(line.rstrip("\n"))
    flush()
    return header, table


def get_samples(header):
    lines = header.get("!Sample_geo_accession", [])
    if not lines:
        return []
    return [c.strip().strip('"') for c in lines[0].split("\t")]


def build_meta(header, samples):
    recs = []
    for i, sid in enumerate(samples):
        rec = {"sample_id": sid}
        for field, col in [
            ("!Sample_title", "Sample_title"),
            ("!Sample_source_name_ch1", "Sample_source_name_ch1"),
        ]:
            rec[col] = ""
            if field in header:
                for ln in header[field]:
                    parts = ln.split("\t")
                    if i < len(parts):
                        v = parts[i].strip().strip('"')
                        if v:
                            rec[col] = v
                            break
        chars = []
        if "!Sample_characteristics_ch1" in header:
            for ln in header["!Sample_characteristics_ch1"]:
                parts = ln.split("\t")
                if i < len(parts):
                    c = parts[i].strip().strip('"')
                    if c:
                        chars.append(c)
        rec["characteristics"] = "; ".join(chars)
        recs.append(rec)
    return pd.DataFrame(recs)


def classify(row):
    text = " ".join(str(v) for v in [
        row.get("Sample_title", ""),
        row.get("Sample_source_name_ch1", ""),
        row.get("characteristics", ""),
    ]).lower()
    for kw in CONTROL_KEYWORDS:
        if kw in text:
            return "control"
    for kw in CASE_KEYWORDS:
        if kw in text:
            return "case"
    return "unknown"


def parse_characteristics(header, samples):
    """返回 (records, key_counter, value_counter)。"""
    rows = header.get("!Sample_characteristics_ch1", [])
    key_counter = Counter()
    value_counter = defaultdict(Counter)
    records = []
    for ln in rows:
        parts = ln.split("\t")
        for i in range(min(len(parts), len(samples))):
            cell = parts[i].strip().strip('"')
            if not cell:
                continue
            if ":" in cell:
                key, _, value = cell.partition(":")
                key = key.strip()
                value = value.strip()
            else:
                key = "(no key)"
                value = cell.strip()
            records.append((samples[i], key, value))
            key_counter[key] += 1
            value_counter[key][value] += 1
    return records, key_counter, value_counter


def parse_table(table_lines):
    header = table_lines[0].split("\t")
    header = [h.strip().strip('"') for h in header]
    probe_col = header[0]
    rows = []
    for line in table_lines[1:]:
        if not line.strip():
            continue
        cells = [c.strip().strip('"') for c in line.split("\t")]
        if len(cells) != len(header):
            cells = cells[:len(header)] + [""] * max(0, len(header) - len(cells))
        rows.append(cells)
    df = pd.DataFrame(rows, columns=header)
    df[probe_col] = df[probe_col].astype(str)
    df = df.set_index(probe_col)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def scale_summary(df):
    vals = df.values.ravel()
    vals = vals[~np.isnan(vals)]
    return {
        "min": float(np.min(vals)),
        "median": float(np.median(vals)),
        "max": float(np.max(vals)),
        "p1": float(np.percentile(vals, 1)),
        "p99": float(np.percentile(vals, 99)),
    }


def fmt_scale(d):
    return (f"min={d['min']:.4f}  median={d['median']:.4f}  max={d['max']:.4f}  "
            f"1%={d['p1']:.4f}  99%={d['p99']:.4f}")


def main():
    log = ["# P1 解析日志", ""]

    for gse, fname in DATASETS.items():
        print(f"\n===== {gse} =====")
        path = RAW / fname
        if not path.exists():
            print(f"[missing] {path}")
            log.append(f"## {gse}\n\n- 状态：源文件缺失，未解析\n")
            continue

        lines = read_lines(path)
        header, table = parse_series_matrix(lines)
        samples = get_samples(header)
        print(f"样本数：{len(samples)}")

        meta = build_meta(header, samples)
        meta["group"] = meta.apply(classify, axis=1)
        meta["split"] = "discovery" if gse == "GSE66360" else "validation"

        # characteristics 键值统计（原样打印）
        records, key_counter, value_counter = parse_characteristics(header, samples)
        print(f"\n--- {gse} characteristics 键值频次 ---")
        for key, cnt in key_counter.most_common():
            print(f"键 [{key}] 出现 {cnt} 次")
            for val, vc in value_counter[key].most_common():
                print(f"    取值 [{val}] -> {vc} 次")

        # 尺度探测
        expr = parse_table(table)
        sc = scale_summary(expr)
        print(f"\n--- {gse} 尺度探测（原始数值矩阵）---")
        print(fmt_scale(sc))
        is_log2 = sc["max"] <= 100
        scale_verdict = "已 log2 转换" if is_log2 else "线性尺度（未 log）"
        print(f"判定：max={sc['max']:.4f} {'<=100' if is_log2 else '>100'} -> {scale_verdict}")

        # 样本来源
        src_values = meta["Sample_source_name_ch1"].astype(str).str.strip()
        src_unique = [v for v in src_values.unique() if v and v != ""]
        print(f"\n--- {gse} 样本来源（!Sample_source_name_ch1）---")
        for v in src_unique:
            print(v)

        title_values = meta["Sample_title"].astype(str).str.strip()
        title_unique = [v for v in title_values.unique() if v and v != ""]
        print(f"\n--- {gse} 样本标题（!Sample_title）---")
        for v in title_unique:
            print(v)

        # 检测 title 中的 DISCOVERY / VALIDATION 子集标注（仅报告，不改变按数据集划分的 split）
        arm = meta["Sample_title"].astype(str).str.extract(r"(DISCOVERY|VALIDATION)", expand=False)
        meta["title_arm"] = arm.fillna("").str.upper()
        arm_counts = meta["title_arm"].value_counts()
        print(f"\n--- {gse} title 中 DISCOVERY/VALIDATION 标注 ---")
        if (meta["title_arm"] != "").any():
            for a, c in arm_counts.items():
                if a:
                    print(f"{a}: {c}")
        else:
            print("未找到")

        # role 映射（P1b）：GSE66360 内部按 title_arm 分 train/internal_val；GSE48060 全 external
        if gse == "GSE66360":
            meta["role"] = meta["title_arm"].map(
                {"DISCOVERY": "train", "VALIDATION": "internal_val"}
            ).fillna("train")
        else:
            meta["role"] = "external"
        print(f"\n--- {gse} role 统计 ---")
        for r, c in meta["role"].value_counts().items():
            print(f"{r}: {c}")

        # group 统计
        group_counts = meta["group"].value_counts()
        print(f"\n--- {gse} group 统计 ---")
        for k, v in group_counts.items():
            print(f"{k}: {v}")
        unknown = meta[meta["group"] == "unknown"]
        if len(unknown):
            print(f"无法判断的样本（{len(unknown)} 个）：")
            print(unknown["sample_id"].tolist())

        # 保存
        out_expr = INTERIM / f"{gse}_expr.parquet"
        out_meta = INTERIM / f"{gse}_meta.csv"
        expr.to_parquet(out_expr)
        meta.to_csv(out_meta, index=False, encoding="utf-8-sig")
        print(f"\n[写] {out_expr.name}  ({expr.shape[0]} 探针 x {expr.shape[1]} 样本)")
        print(f"[写] {out_meta.name}")

        # 写入 log
        log.append(f"## {gse}\n")
        log.append(f"### 数据集概况\n")
        log.append(f"- 样本数：{len(samples)}")
        log.append(f"- 样本来源（!Sample_source_name_ch1 唯一值）：{src_unique if src_unique else '无'}")
        log.append(f"- 样本标题（!Sample_title 唯一值）：{title_unique if title_unique else '无'}")
        log.append("")
        log.append("### 分组映射\n")
        log.append("- 规则：case <- `acute myocardial infarction` / `stemi` / `ami`；control <- `healthy` / `normal` / `control`")
        for k, v in group_counts.items():
            log.append(f"- {k}: {v}")
        if len(unknown):
            log.append(f"- 无法判断样本（未猜，单独列出）：{unknown['sample_id'].tolist()}")
        else:
            log.append("- 无法判断样本：无")
        log.append("")
        log.append("### 尺度探测\n")
        log.append(f"- {fmt_scale(sc)}")
        log.append(f"- 判定：{scale_verdict}（{sc['max']:.4f} {'<=100' if is_log2 else '>100'}）")
        log.append(f"- 说明：本阶段不做任何变换，仅记录。")
        log.append("")
        log.append("### characteristics 键值频次\n")
        for key, cnt in key_counter.most_common():
            log.append(f"- 键 `{key}`（{cnt} 次）")
            for val, vc in value_counter[key].most_common():
                log.append(f"  - `{val}` -> {vc}")
        log.append("")
        log.append("### discovery/validation 检测\n")
        log.append("- characteristics 中未找到 discovery/validation 字段。")
        if (meta["title_arm"] != "").any():
            log.append(f"- **!Sample_title 中含官方 DISCOVERY/VALIDATION 标注**：{dict(arm_counts.dropna())}")
            log.append("- 按既定约定，split 列仍按数据集划分：GSE66360=discovery，GSE48060=validation。")
        else:
            log.append("- split 按数据集划分：GSE66360=discovery，GSE48060=validation（外部验证集，不参与调参）。")
        log.append("")

    # 探针合并规则（预留，本阶段不执行）
    log.append("## 探针合并规则（预留）\n")
    log.append("- 本阶段（P1）不执行多探针合并。")
    log.append("- 后续 03_preprocess.py 将实现：多探针映射同一基因时取 **median**。")
    log.append("- 03 脚本已预留 `merge_probes()` 实现位置（P2 提供）。")
    log.append("")

    # 异常字段
    log.append("## 异常字段说明\n")
    log.append("- 本阶段未做过滤，矩阵中存在探针缺失/非数值时以 NaN 记录，待 P2 处理。")
    log.append("")

    # 写 parse_log.md
    (LOGS / "parse_log.md").write_text("\n".join(log), encoding="utf-8")

    print(f"\n[02_parse_geo] 已写入 {LOGS / 'parse_log.md'}")


if __name__ == "__main__":
    main()
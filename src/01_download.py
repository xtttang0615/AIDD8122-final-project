# -*- coding: utf-8 -*-
"""下载 GEO series matrix 与 GPL570 注释到 data/raw/，并生成 MD5 校验。

输入：无（网络下载）
输出：data/raw/*.gz、data/raw/checksums.md5
运行：C:/Anaconda/Anaconda3/python.exe src/01_download.py
"""
from pathlib import Path
import hashlib
import gzip
import sys
import time

import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

FILES = {
    "GSE66360_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66360/matrix/GSE66360_series_matrix.txt.gz",
    "GSE48060_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48060/matrix/GSE48060_series_matrix.txt.gz",
    "GPL570.annot.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
}

MIN_SIZE = 1_000_000  # 1 MB
TIMEOUT = 60          # 连接/读取超时（秒）
RETRIES = 3
MAX_DURATION = 600    # 单文件最长 10 分钟


def file_md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_ok(path: Path) -> bool:
    """读完整文件，确认 gzip 可解压。"""
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except Exception:
        return False


def download(url: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > MIN_SIZE:
        print(f"[skip] {dest.name} ({dest.stat().st_size} bytes)")
        return True

    for attempt in range(1, RETRIES + 1):
        print(f"[get] {dest.name}  attempt {attempt}/{RETRIES}")
        start = time.monotonic()
        tmp = dest.with_name(dest.name + ".part")
        try:
            with requests.get(url, stream=True, timeout=TIMEOUT) as r:
                r.raise_for_status()
                total = int(r.headers.get("Content-Length") or 0)
                done = 0
                with tmp.open("wb") as f:
                    for chunk in r.iter_content(chunk_size=256 * 1024):
                        if time.monotonic() - start > MAX_DURATION:
                            print(f"\n[stop] {dest.name} 超过 10 分钟无进展，请手动浏览器下载并保存到 data/raw/{dest.name}")
                            return False
                        if not chunk:
                            continue
                        f.write(chunk)
                        done += len(chunk)
                        if total:
                            sys.stdout.write(f"\r  {done * 100 // total}% ({done}/{total} bytes)")
                            sys.stdout.flush()
                print()
            if not gzip_ok(tmp):
                raise ValueError("gzip 解压校验失败")
            tmp.replace(dest)
            print(f"[ok] {dest.name} -> {dest.stat().st_size} bytes")
            return True
        except Exception as e:
            print(f"[retry] {dest.name} attempt {attempt} 失败: {e}")
            if tmp.exists():
                tmp.unlink()
            if attempt == RETRIES:
                print(f"[fail] {dest.name} 重试 {RETRIES} 次仍失败，请手动浏览器下载并保存到 data/raw/{dest.name}")
                return False
            time.sleep(2)
    return False


def main():
    results = {}
    for name, url in FILES.items():
        results[name] = download(url, RAW / name)

    md5_lines = []
    for name in FILES:
        p = RAW / name
        if p.exists():
            md5_lines.append(f"{file_md5(p)}  {name}")
    (RAW / "checksums.md5").write_text("\n".join(md5_lines) + "\n", encoding="utf-8")

    ok = sum(1 for v in results.values() if v)
    print(f"[01_download] {ok}/{len(FILES)} 文件就绪；checksums.md5 已写入")
    for name, v in results.items():
        if not v:
            print(f"[01_download] 需手工处理: {name}")


if __name__ == "__main__":
    main()
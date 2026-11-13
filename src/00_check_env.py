# -*- coding: utf-8 -*-
"""
输入文件：无
输出文件：logs/env_report.md、setup/missing_packages.txt
运行命令：python src/00_check_env.py
说明：自检 Python 侧系统信息与科学计算包版本，不安装任何包。
"""
from pathlib import Path
import os
import sys
import platform
import ctypes
import shutil
import subprocess
from datetime import datetime
from importlib import metadata as md

ROOT = Path(__file__).resolve().parents[1]
LOG_DIR = ROOT / "logs"
SETUP_DIR = ROOT / "setup"
LOG_DIR.mkdir(parents=True, exist_ok=True)
SETUP_DIR.mkdir(parents=True, exist_ok=True)
REPORT = LOG_DIR / "env_report.md"
MISSING = SETUP_DIR / "missing_packages.txt"

# 需要检查的包（label 用于显示，dist 是 PyPI 发行名）
PACKAGES = [
    ("numpy", "numpy"),
    ("pandas", "pandas"),
    ("scipy", "scipy"),
    ("sklearn", "scikit-learn"),
    ("matplotlib", "matplotlib"),
    ("seaborn", "seaborn"),
    ("rdkit", "rdkit"),
    ("xgboost", "xgboost"),
    ("statsmodels", "statsmodels"),
    ("pyarrow", "pyarrow"),
    ("requests", "requests"),
]


def cpu_cores():
    return os.cpu_count()


def total_memory_gb():
    """Windows 下用 GlobalMemoryStatusEx 读取总物理内存（GB）。"""
    try:
        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_uint32),
                ("dwMemoryLoad", ctypes.c_uint32),
                ("ullTotalPhys", ctypes.c_uint64),
                ("ullAvailPhys", ctypes.c_uint64),
                ("ullTotalPageFile", ctypes.c_uint64),
                ("ullAvailPageFile", ctypes.c_uint64),
                ("ullTotalVirtual", ctypes.c_uint64),
                ("ullAvailVirtual", ctypes.c_uint64),
                ("ullAvailExtendedVirtual", ctypes.c_uint64),
            ]
        m = MEMORYSTATUSEX()
        m.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
        return round(m.ullTotalPhys / (1024 ** 3), 1)
    except Exception:
        return None


def gpu_summary():
    """检测 nvidia-smi 是否可用，可用则读取 GPU 名称与显存。"""
    nv = shutil.which("nvidia-smi")
    if not nv:
        return "not found (nvidia-smi 不在 PATH 中)"
    try:
        out = subprocess.run(
            [nv, "--query-gpu=name,memory.total", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=15,
        )
        if out.returncode == 0 and out.stdout.strip():
            return out.stdout.strip().replace("\r", "")
        return "nvidia-smi found but query failed"
    except Exception:
        return "nvidia-smi found but not callable"


def pkg_version(dist):
    try:
        return md.version(dist)
    except Exception:
        return "MISSING"


def main():
    lines = []
    lines.append("# Environment Report (Python + R)")
    lines.append("")
    lines.append(f"> generated at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    lines.append("## 1. System")
    lines.append("")
    lines.append(f"- OS: {platform.system()} {platform.release()} ({platform.version()})")
    lines.append(f"- Machine: {platform.machine()}")
    lines.append(f"- CPU cores: {cpu_cores()}")
    lines.append(f"- Total memory (GB): {total_memory_gb()}")
    lines.append(f"- NVIDIA GPU: {gpu_summary()}")
    lines.append("")
    lines.append("## 2. Python")
    lines.append("")
    lines.append(f"- version: {sys.version.split()[0]}")
    lines.append(f"- executable: {sys.executable}")
    lines.append("")
    lines.append("## 3. Python packages")
    lines.append("")
    lines.append("| package | version |")
    lines.append("|---|---|")
    missing = []
    for label, dist in PACKAGES:
        v = pkg_version(dist)
        lines.append(f"| {label} | {v} |")
        if v == "MISSING":
            missing.append(dist)
    lines.append("")

    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    if missing:
        MISSING.write_text("Python packages missing:\n" + "\n".join(missing) + "\n", encoding="utf-8")
    else:
        MISSING.write_text("Python packages missing: none\n", encoding="utf-8")

    print(f"[00_check_env.py] report written to {REPORT}")
    print(f"[00_check_env.py] missing packages: {missing if missing else 'none'}")


if __name__ == "__main__":
    main()
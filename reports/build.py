import base64, os, re
from pathlib import Path

HERE = Path(__file__).resolve().parent          # .../reports
ROOT = HERE.parent                               # project / release root
BASE = [ROOT / "results" / "figures", HERE, ROOT]

parts = [HERE / f"part{i}.html" for i in range(1, 5)]
html = "\n".join(open(p, encoding="utf-8").read() for p in parts)
html = html.replace("</body>\n</html>", "")
html = html.replace("</body>", "").replace("</html>", "")


def repl(m):
    name = m.group(1)
    for b in BASE:
        p = Path(b) / name
        if p.exists():
            data = base64.b64encode(open(p, "rb").read()).decode()
            return f"data:image/png;base64,{data}"
    raise FileNotFoundError(name)


html = re.sub(r"__IMG__(.+?)__", repl, html)
html += "\n</body>\n</html>\n"
out = HERE / "report_full.html"
out.write_text(html, encoding="utf-8")
print("wrote", out, "size MB", round(len(html) / 1e6, 2))

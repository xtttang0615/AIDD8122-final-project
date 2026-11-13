# 文档检查记录（DOCX / PDF）

## 环境限制（必须如实记录）

- 本机**没有安装 Microsoft Word**，也没有 LibreOffice（`soffice` 不在 PATH）。
- 因此**无法**用 Word COM 或 LibreOffice 更新 TOC 域、无法重排分页、无法逐页渲染截图检查。
- 未安装 `python-docx`（未擅自安装），未做 DOCX 内部样式后处理。

## 已完成的处理

| 项目 | 处理 |
|---|---|
| 正文来源 | 唯一来源为 `reports/part1–part4.html` → `reports/build.py` → `report_full.html`（图片全部 base64 内联） |
| DOCX 生成 | `pandoc report_full.html --toc --toc-depth=2`，生成 `reports/final/AIDD8122_Final_Submission.docx` |
| 目录 | pandoc 写入真实 `TOC \o "1-2"` 域（已在 `word/document.xml` 中确认），**首次打开需 F9 更新** |
| PDF 生成 | Edge headless `--print-to-pdf`，无浏览器页眉页脚、无 `file:///` 页脚 |
| 图片 | 全部内联为 data URI，DOCX 内嵌图片，无外部链接 |

## 未完成（因环境限制）

- TOC 页码未更新（域存在但为未计算状态）。
- 未逐页渲染检查：跨页表格、图表分页、图注是否与图分离、空白页、字号可读性。
- 未核对图表编号与交叉引用的最终一致性（仅在 HTML 层面确认编号连续）。

## 交付时的用户操作

1. 用 Word 打开 `reports/final/AIDD8122_Final_Submission.docx`。
2. `Ctrl+A` → `F9` → 选择“更新整个目录”。
3. 检查目录页码与正文一致后另存；如需 PDF，用 Word「另存为 PDF」以保留分页。
4. 封面的 GitHub 链接仍为占位符，需替换。

## 内容一致性核查（本轮）

- 表7 第 1–10 名：已逐行核对，源于 eports/table7_check.csv（mean-score 与 PDE4B-only 两个列表，含 ChEMBL ID、三靶点得分、prob_mean、排名）。
- 重叠数值：eports/top_overlap_check.csv。比较对象明确为：
  mean_score vs mean_percentile = 11/20（30/50）；mean_percentile vs PDE4B_only = 5/20（15/50）；
  mean_score vs PDE4B_only = 12/20（31/50）；mean_percentile vs Borda = 20/20（49/50）。
  固定取前 N，不扩展同分；分母 1,823。
- 附录 B：条目 1–13 标注为历史问题及已实施修正；新增 14–18（身份连接、训练负例重叠、统一骨架分折、no-op 消融、PDK4 重复折）。
- DOCX 内容已在 word/document.xml 中验证含本轮修改（1,282 / Entries 1 / unconfirmed PDE4B activity / 12 of 20 / 5 of 20 / 1,823）。
- 仍未执行：目录域页码计算、逐页视觉检查（无 Word / LibreOffice）。

## 附录 A 与参考文献核查（src/31_verify_appendix_refs.py）

- **附录 A 五折基因列表**：由 logs/model_log.md 重新解析并与 part4.html 附录 A 表格逐折比对，**5/5 折完全一致**（ALL_MATCH True）。跨折频率由同一份列表计算，与附录中的频率表同源。
- **引用编号**：正文引用编号集合与文献表编号集合完全相同，共 24 条，无缺号、无重复、最大编号 24。
- **条目完整性**：24 条中 21 条含 DOI 与年份；[19] JMLR（该刊无 DOI）、[21] RDKit（软件，给 URL）、[23] SciPy 会议（该条按会议格式引用）无 DOI，属格式正常，非缺失。
- **未核**：文献与原文的对应内容（作者、卷期、结论）未逐条比对原文，属外部事实核对，不在本仓库可验证范围。

## 参考文献真实性核验（Crossref API，src/33_verify_references.py）

方法：抽取决议条目中的 DOI，逐条调用 pi.crossref.org/works/{doi}，比对返回标题与本文引用标题的词元重合度。

| 结果 | 条目 |
|---|---|
| OK（DOI 可解析且标题匹配） | [1][3][4][5][6][7][8][9][10][11][12][13][14][15][16][17][18][20][22][24] 共 20 条 |
| **已修正** | [2] DOI 原为 10.1016/j.yjmcc.2014.04.009（解析为 "Elucidating ERK2 function in the heart"，无关文献）；经 Crossref 标题检索更正为 10.1016/j.yjmcc.2014.04.017，标题、期刊（J Mol Cell Cardiol 2014;74:13-21）均匹配 |
| 无 DOI（格式正常，非缺失） | [19] JMLR（该刊不赋 DOI）；[21] RDKit 软件（给官网 URL）；[23] SciPy 会议论文（会议格式） |

结果表：eports/reference_check.csv。

**未核**：各文献结论是否支持正文具体论断（如 [12] 的 PDE4B 心脏保护方向、[15] 的 ibudilast 优先抑制 PDE4），仅核到标题/期刊/年份层面。

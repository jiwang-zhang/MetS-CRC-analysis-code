#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE10950 独立队列验证：7 hub 基因在 CRC 中的差异表达可重复性。

- 输入：
    GSE10950.csv  (GeneName + GSMxxx_con / GSMxxx_tre 交替列；con=正常, tre=肿瘤)
    主分析 co-DEGs_CRC_MetS.csv (含 CRC_log2FC 方向，用于方向一致性对比)
- 方法：
    对每个 hub 基因，肿瘤(_tre) vs 正常(_con)：
      * log2FC = mean(tre) - mean(con)   (矩阵按 GEO GPL6104 惯例为 log2 尺度)
      * Welch t-test (scipy ttest_ind equal_var=False)
      * Mann-Whitney U (非参)
      * 7 基因 BH 校正 (FDR)
    Direction consistency：sign(GSE10950 log2FC) == sign(主分析 CRC_log2FC)
- 输出：
    GSE10950_hub_validation.csv
    GSE10950_hub_validation.png  (log2FC 柱状图，颜色标注显著/一致)
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"D:\BaiduSyncdisk\daixie-colorectal\第二部分DEG"
OUT  = os.path.join(BASE, "27.GSE10950_validation")
os.makedirs(OUT, exist_ok=True)

GSE_FILE  = os.path.join(BASE, "25.diagnostic_model", "GSE10950.csv")
CO_FILE   = os.path.join(BASE, "09.两个疾病分别DEG", "co_DEGs_CRC_MetS.csv")
HUBS = ["TNFRSF21", "LPCAT4", "PLAU", "GPAT3", "SGK1", "TMCC3", "PPARD"]

# ---------- 1. 载入 GSE10950 ----------
df = pd.read_csv(GSE_FILE)
df = df.set_index("GeneName")
con_cols = [c for c in df.columns if c.endswith("_con")]
tre_cols = [c for c in df.columns if c.endswith("_tre")]
assert len(con_cols) == len(tre_cols), f"con/tre 数量不匹配: {len(con_cols)} vs {len(tre_cols)}"
print(f"[info] GSE10950: {df.shape[0]} 基因, {len(con_cols)} 正常 / {len(tre_cols)} 肿瘤")

# ---------- 2. 主分析 CRC 方向 ----------
co = pd.read_csv(CO_FILE)
co_dir = dict(zip(co["Gene"], co["CRC_log2FC"]))

# ---------- 3. 逐基因检验 ----------
rows = []
missing = []
for g in HUBS:
    if g not in df.index:
        missing.append(g); continue
    y_tre = df.loc[g, tre_cols].to_numpy(dtype=float)
    y_con = df.loc[g, con_cols].to_numpy(dtype=float)
    mean_tre, mean_con = y_tre.mean(), y_con.mean()
    log2fc = mean_tre - mean_con
    t, p_t = stats.ttest_ind(y_tre, y_con, equal_var=False)
    u, p_u = stats.mannwhitneyu(y_tre, y_con, alternative="two-sided")
    main_dir = co_dir.get(g, np.nan)
    consistent = (np.sign(log2fc) == np.sign(main_dir)) if (not np.isnan(main_dir) and main_dir != 0) else np.nan
    rows.append({
        "Gene": g,
        "mean_normal": round(mean_con, 4),
        "mean_tumor":  round(mean_tre, 4),
        "log2FC_tumor_vs_normal": round(log2fc, 4),
        "Welch_t": round(t, 4),
        "p_Welch": p_t,
        "p_MannWhitney": p_u,
        "main_CRC_log2FC": round(main_dir, 4) if not np.isnan(main_dir) else np.nan,
        "direction_consistent_with_main": consistent,
    })

res = pd.DataFrame(rows)
# BH 校正（基于 Welch p）
res["p_Welch_FDR"] = stats.false_discovery_control(res["p_Welch"].values, method="bh")
res["p_MannWhitney_FDR"] = stats.false_discovery_control(res["p_MannWhitney"].values, method="bh")
res["significant_FDR<0.05"] = res["p_Welch_FDR"] < 0.05

# 排序：保持 hub 顺序
out_csv = os.path.join(OUT, "GSE10950_hub_validation.csv")
res.to_csv(out_csv, index=False)
print(f"\n[result] 写入 {out_csv}\n")
print(res[["Gene","log2FC_tumor_vs_normal","p_Welch","p_Welch_FDR",
           "main_CRC_log2FC","direction_consistent_with_main","significant_FDR<0.05"]]
      .to_string(index=False))

if missing:
    print(f"\n[warn] 未在 GSE10950 中找到: {missing}")

# ---------- 4. 绘图 ----------
fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(res))
colors = []
for _, r in res.iterrows():
    if r["significant_FDR<0.05"] and r["direction_consistent_with_main"]:
        colors.append("#2ca02c")   # 绿：显著且方向一致
    elif r["significant_FDR<0.05"]:
        colors.append("#d62728")   # 红：显著但方向不一致（需警惕）
    else:
        colors.append("#7f7f7f")   # 灰：不显著
bars = ax.bar(x, res["log2FC_tumor_vs_normal"], color=colors, edgecolor="black", linewidth=0.6)
ax.axhline(0, color="black", linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels(res["Gene"], rotation=30, ha="right", fontsize=10)
ax.set_ylabel("log2FC (tumor vs normal)", fontsize=11)
ax.set_title("GSE10950 independent validation:\n7 ML-identified hub genes in CRC", fontsize=12)
for xi, (_, r) in zip(x, res.iterrows()):
    sig = "*" if bool(r["significant_FDR<0.05"]) else ""
    ax.text(xi, r["log2FC_tumor_vs_normal"] + (0.05 if r["log2FC_tumor_vs_normal"] >= 0 else -0.05),
            f"{r['log2FC_tumor_vs_normal']:+.2f}{sig}", ha="center",
            va="bottom" if r["log2FC_tumor_vs_normal"] >= 0 else "top", fontsize=9)
from matplotlib.patches import Patch
legend = [Patch(facecolor="#2ca02c", label="Significant & consistent"),
          Patch(facecolor="#d62728", label="Significant but discordant"),
          Patch(facecolor="#7f7f7f", label="Not significant")]
ax.legend(handles=legend, fontsize=8, loc="best")
plt.tight_layout()
out_png = os.path.join(OUT, "GSE10950_hub_validation.png")
fig.savefig(out_png, dpi=300)
print(f"\n[figure] 写入 {out_png}")

# ---------- 5. 一致性汇总 ----------
n_sig = int(res["significant_FDR<0.05"].sum())
n_cons = int((res["significant_FDR<0.05"] & res["direction_consistent_with_main"]).sum())
print(f"\n[summary] 显著(FDR<0.05): {n_sig}/{len(res)}; 其中方向一致: {n_cons}")

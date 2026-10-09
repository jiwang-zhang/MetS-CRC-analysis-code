"""
MetS 血液 RNA-seq 差异分析
  - GSE200744: 原始计数 -> PyDESeq2 (MetS 25 vs Normal 33)
  - GSE145412: 归一化表达 -> limma 思路 empirical-Bayes moderated t (MetS 16 vs Healthy 16, 敏感性)
阈值: |log2FC| > 0.5 且 padj < 0.05 (与模板一致)
"""
import os, re, gzip, csv
import numpy as np
import pandas as pd

BASE = r"D:\daixie-colorectal\19.databases"
OUT_DIR = BASE

# ----------------------------------------------------------------------
# 工具: limma 风格 empirical-Bayes moderated t-test (两组的简化实现)
# ----------------------------------------------------------------------
def limma_modt(expr, grp):
    """expr: genes x samples (numpy, 已 log2 变换); grp: 0/1 数组(1=case)
    返回 per-gene: log2FC, mod_t, pvalue, df_total"""
    expr = np.asarray(expr, float)
    g1 = grp == 1
    g0 = grp == 0
    n1, n0 = g1.sum(), g0.sum()
    m1 = expr[:, g1].mean(1)
    m0 = expr[:, g0].mean(1)
    # 各基因组内残差方差(等方差假设)
    ss1 = ((expr[:, g1] - m1[:, None]) ** 2).sum(1)
    ss0 = ((expr[:, g0] - m0[:, None]) ** 2).sum(1)
    df_res = n1 + n0 - 2
    s2 = (ss1 + ss0) / df_res
    s2 = np.where(s2 <= 0, 1e-8, s2)
    fc = m1 - m0
    se = np.sqrt(s2 * (1.0 / n1 + 1.0 / n0))
    t_raw = fc / se
    # empirical Bayes 收缩 (limma squeezeVar, method-of-moments)
    # 先验方差 s20 与先验自由度 d0
    s2_med = np.median(s2)
    # 用经验贝叶斯估计 d0, s20 (近似 limma)
    # 先估计 s20 = median(s2); d0 由方差的变异估计
    var_s2 = np.var(s2)
    d0 = (2 * s2_med**2) / max(var_s2, 1e-8) if var_s2 > 0 else 100.0
    d0 = min(max(d0, 1.0), 1000.0)
    s20 = s2_med
    s2_post = (d0 * s20 + df_res * s2) / (d0 + df_res)
    mod_t = fc / np.sqrt(s2_post * (1.0 / n1 + 1.0 / n0))
    df_tot = d0 + df_res
    from scipy import stats as st
    pval = 2 * st.t.sf(np.abs(mod_t), df_tot)
    return fc, mod_t, pval, df_tot, (n1, n0)

def bh(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    # NaN 视为不显著(置为1), 避免 minimum.accumulate 把整列传染成 NaN
    p = np.where(np.isnan(p), 1.0, p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(1, n + 1))
    # 从后往前取极小值，保证单调不增
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out

# ----------------------------------------------------------------------
# 1) GSE200744 (DESeq2)
# ----------------------------------------------------------------------
def deg_gse200744():
    print("\n" + "=" * 70)
    print("GSE200744  DESeq2: MetS vs Normal")
    f = os.path.join(BASE, "GSE200744", "GSE200744_RNAseq_count_in_Mets.txt.gz")
    ens_col = "Ensembl Gene ID"
    sym_col = "Gene Symbol"
    df = pd.read_csv(f, sep="\t", compression="gzip")
    df = df.set_index(ens_col)                 # 以 Ensembl 为基因标识
    df = df[~df.index.duplicated(keep="first")]
    sample_cols = [c for c in df.columns if re.match(r"^(mets|premets|normal)_\d+$", c)]
    def grp_of(c):
        if c.startswith("mets_"): return "MetS"
        if c.startswith("premets_"): return "PreMetS"
        return "Normal"
    sample_cols2 = [c for c in sample_cols if grp_of(c) in ("MetS", "Normal")]
    counts = df[sample_cols2].astype(int)
    # 低表达过滤
    counts = counts[counts.sum(axis=1) >= 10]
    meta = pd.DataFrame({"condition": [grp_of(c) for c in counts.columns]},
                        index=counts.columns)
    meta["condition"] = pd.Categorical(meta["condition"], categories=["Normal", "MetS"])
    print(f"  样本: MetS={ (meta.condition=='MetS').sum() }, Normal={ (meta.condition=='Normal').sum() }")
    print(f"  基因数(过滤后)={counts.shape[0]}")
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    # pydeseq2(anndata) 要求 counts 为 样本×基因 -> 转置
    dds = DeseqDataSet(counts=counts.T, metadata=meta, design="~condition", refit_cooks=True)
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=["condition", "MetS", "Normal"])
    stat_res.summary()
    res = stat_res.results_df.copy()
    res.index.name = "Ensembl"
    res = res.reset_index()
    res = res.rename(columns={"log2FoldChange": "log2FC",
                              "padj": "adj.P.Val", "pvalue": "P.Value"})
    sym_map = dict(zip(df.index, df[sym_col]))
    res["Gene Symbol"] = [str(sym_map.get(e, "")) for e in res["Ensembl"]]
    res = res[["Ensembl", "Gene Symbol", "baseMean", "log2FC", "lfcSE", "stat", "P.Value", "adj.P.Val"]]
    sig = res[(res["adj.P.Val"] < 0.05) & (res["log2FC"].abs() > 0.5)]
    up = sig[sig["log2FC"] > 0]; down = sig[sig["log2FC"] < 0]
    print(f"  DEG (|log2FC|>0.5 & padj<0.05): {len(sig)}  (上调 {len(up)}, 下调 {len(down)})")
    out = os.path.join(OUT_DIR, "GSE200744", "DEG_MetS_vs_Normal.csv")
    res.to_csv(out, index=False)
    sig.to_csv(os.path.join(OUT_DIR, "GSE200744", "DEG_MetS_vs_Normal.sig.csv"), index=False)
    print(f"  写出: {out}")
    return res, sig

# ----------------------------------------------------------------------
# 2) GSE145412 (limma 思路)
# ----------------------------------------------------------------------
def map_gse145412_samples():
    """series_matrix title 解析 sample_X -> (health, bmi)"""
    p = os.path.join(BASE, "GSE145412", "GSE145412_series_matrix.txt.gz")
    acc = title = None
    with gzip.open(p, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                acc = line.rstrip("\n").split("\t")[1:]
            elif line.startswith("!Sample_title"):
                title = line.rstrip("\n").split("\t")[1:]
    assert acc and title, "解析失败"
    mapping = {}
    for a, t in zip(acc, title):
        a = a.strip().strip('"'); t = t.strip().strip('"')
        m = re.search(r"sample_(\d+)", t)
        if not m:
            continue
        sidx = m.group(1)
        toks = [x.strip() for x in t.split(",")]
        health = "MetS" if "D" in toks else ("Healthy" if "H" in toks else "NA")
        bmi = "Obese" if "O" in toks else ("Lean" if "L" in toks else "NA")
        mapping[sidx] = (health, bmi)
    return mapping

def deg_gse145412():
    print("\n" + "=" * 70)
    print("GSE145412  limma式 moderated t: MetS vs Healthy (敏感性)")
    f = os.path.join(BASE, "GSE145412", "GSE145412_rnaseq.txt.gz")
    df = pd.read_csv(f, sep="\t", compression="gzip")
    samp_cols = [c for c in df.columns if re.match(r"^sample_\d+$", c)]
    mapping = map_gse145412_samples()
    # 仅保留能映射的样本
    valid = [c for c in samp_cols if c.split("_")[1] in mapping]
    health = [mapping[c.split("_")[1]][0] for c in valid]
    print("  各组样本数:", {k: health.count(k) for k in set(health)})
    expr = df[valid].copy().astype(float)
    expr_log = np.log2(expr.values + 1.0)
    # 低表达过滤: 仅在 >=50% 样本有表达(原始值>0)的基因上检验, 降低多重检验负担
    present = (expr.values > 0).sum(1)
    keep = present >= 0.5 * len(valid)
    expr_log = expr_log[keep]
    print(f"  基因数(低表达过滤后)={expr_log.shape[0]} (过滤前={keep.shape[0]})")
    grp = np.array([1 if h == "MetS" else 0 for h in health])
    fc, modt, pval, dft, (n1, n0) = limma_modt(expr_log, grp)
    padj = bh(pval)
    out = pd.DataFrame({
        "geneID": df["geneID"].values[keep],
        "External_gene_name": df["External_gene_name"].values[keep],
        "EntrezID": df["EntrezID"].values[keep],
        "log2FC": fc, "modt": modt, "P.Value": pval, "adj.P.Val": padj,
    })
    sig = out[(out["adj.P.Val"] < 0.05) & (out["log2FC"].abs() > 0.5)]
    sugg = out[(out["P.Value"] < 0.05) & (out["log2FC"].abs() > 0.5)]   # 提示性(未校正, 供交集)
    up = sig[sig["log2FC"] > 0]; down = sig[sig["log2FC"] < 0]
    print(f"  DEG (BH adj.P.Val<0.05 & |log2FC|>0.5): {len(sig)}  (上调 {len(up)}, 下调 {len(down)})")
    print(f"  提示性 (raw P<0.05 & |log2FC|>0.5): {len(sugg)}")
    out.to_csv(os.path.join(OUT_DIR, "GSE145412", "DEG_MetS_vs_Healthy.csv"), index=False)
    sig.to_csv(os.path.join(OUT_DIR, "GSE145412", "DEG_MetS_vs_Healthy.sig.csv"), index=False)
    sugg.to_csv(os.path.join(OUT_DIR, "GSE145412", "DEG_MetS_vs_Healthy.suggestive.csv"), index=False)
    print(f"  写出: GSE145412/DEG_MetS_vs_Healthy.csv (+ .sig / .suggestive)")
    return out, sig

if __name__ == "__main__":
    r1, s1 = deg_gse200744()
    r2, s2 = deg_gse145412()
    print("\n完成。")

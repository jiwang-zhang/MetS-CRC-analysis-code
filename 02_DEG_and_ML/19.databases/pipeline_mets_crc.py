# -*- coding: utf-8 -*-
"""
MetS x CRC 跨病差异分析管线 (生物信息学阶段, 步骤 3 预演)
==================================================================
本脚本在沙箱(受管 Python venv)中一次跑通 4 个阶段:
  Stage 1 : GSE98895 (GPL6947 Illumina 芯片, PBMC) limma -> MetS vs 对照 DEG
  Stage 2 : MetS 血液共识 = GSE98895 sig ∩ GSE200744 sig (方向一致)
  Stage 3 : CRC 训练集 GPL570 x3 (GSE21510/32323/39582) 各自 limma -> 基因级 IVW meta
  Stage 4 : co-DEGs = MetS 共识 ∩ CRC meta sig (方向一致)

设计要点(对应"教而不替"注释版):
  * 芯片计数在 series_matrix 内嵌;RNA-seq 在 suppl(本脚本不重算, 直接读既定结果)
  * 探针->基因: GPL570.annot.gz 第3列 Gene symbol, 多符号取第一个(///分隔)
  * 同基因多探针: 取均值(mean)聚合
  * 差异: limma 经验贝叶斯 moderatet (自实现, 用方法矩估计先验 d0/s02)
  * 跨研究 meta: 逆方差加权(IVW) 基因级 log2FC, SE=1/sqrt(sum(1/se^2))
  * 阈值统一: |log2FC|>0.5 且 BH adj.P<0.05 (对齐模板)
"""
import os, gzip, io
import numpy as np
import pandas as pd
from scipy import stats

BASE = r"D:/daixie-colorectal/19.databases"
RES  = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)
ANNOT_GPL570  = os.path.join(BASE, "GSE21510", "GPL570.annot.gz")   # CRC 训练集 (Affymetrix HG-U133_Plus_2) 真实 GPL570 注释，置于 CRC 文件夹
ANNOT_GPL6947 = os.path.join(BASE, "GSE98895", "GPL6947.annot.gz")   # GSE98895 真实平台 = Illumina HumanHT-12 V3

# ----------------------------------------------------------------------
# 0. 通用工具
# ----------------------------------------------------------------------
def load_annot(path):
    """解析 GPL570.annot.gz -> {probe_id: gene_symbol}(取第一个///符号)"""
    sym_map = {}
    with gzip.open(path, "rt", errors="ignore") as f:
        header = None
        for line in f:
            if line.startswith("!platform_table_begin"):
                continue
            if line.startswith("#") or line.startswith("!") or line.startswith("^"):
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                i_id = header.index("ID")
                i_sym = header.index("Gene symbol")
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_id, i_sym):
                continue
            pid = parts[i_id]
            sym = parts[i_sym].split("///")[0].strip()
            if sym and sym not in ("", "---", "NA", "NaN"):
                sym_map[pid] = sym
    return sym_map

def parse_meta(path):
    """读取 series_matrix 的 !Sample_* 元数据(只读表头前, 不读大表)"""
    meta = {}
    with gzip.open(path, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_") or line.startswith("!series_"):
                parts = line.rstrip("\n").split("\t")
                meta[parts[0]] = parts[1:]
    return meta

def read_expr(path):
    """流式读取表达表(comment='!' 跳过所有元数据/begin/end 行)"""
    df = pd.read_csv(path, sep="\t", index_col=0, comment="!",
                     compression="gzip", low_memory=False)
    return df

def probes_to_genes(expr, sym_map):
    """index=probe -> 映射为 gene symbol, 同基因多探针取均值"""
    expr = expr.copy()
    expr.index = [sym_map.get(str(p), None) for p in expr.index]
    expr = expr[~expr.index.isna()]
    expr.index = expr.index.astype(str)
    # 数值化
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.groupby(level=0, sort=False).mean()
    return expr

def maybe_log2(expr):
    """GEO 处理数据常见已是 log2; 若数值过大(原始强度)则补 log2。
    Illumina 背景校正会出现负值 -> 加偏移量使其>=1 再 log2, 避免 NaN。"""
    arr = expr.to_numpy()
    arr = arr[np.isfinite(arr)]
    mx, mn = arr.max(), arr.min()
    if mx > 100:   # 原始强度量级 -> log2 变换
        off = (max(0.0, -mn) + 1.0) if mn < 0 else 0.0
        expr = np.log2(expr + off)
        print(f"    [log2] max={mx:.1f} min={mn:.1f} -> log2(+offset {off:.1f})")
    else:
        print(f"    [log2] max={mx:.1f} min={mn:.1f} -> 视为已 log2, 不转换")
    return expr

# ----------------------------------------------------------------------
# 1. limma 经验贝叶斯 moderatet (自实现)
# ----------------------------------------------------------------------
def limma_modt(expr, grp):
    """
    expr: DataFrame, 行=基因(或探针), 列=样本(数值, 已 log2)
    grp : 长度=列数, 0/1 (1=case)
    返回 DataFrame: log2FC, modt, P.Value, SE(修正后), 方向
    """
    grp = np.asarray(grp, dtype=float)
    y = expr.T.values.astype(float)          # n x G
    n, G = y.shape
    # 设计矩阵 [intercept, grp]
    X = np.column_stack([np.ones(n), grp])    # n x 2
    XtX = X.T @ X
    XtX_inv = np.linalg.inv(XtX)
    beta = XtX_inv @ (X.T @ y)               # 2 x G  (beta[1]=grp coef)
    fitted = X @ beta
    resid = y - fitted
    df_res = n - 2
    sigma2 = (resid ** 2).sum(axis=0) / df_res   # G 残差方差
    # 经验贝叶斯先验(方法矩估计, 等 df)
    m = np.mean(sigma2)
    v = np.var(sigma2, ddof=1)
    denom = v - 2.0 * m * m / df_res
    d0 = (2.0 * m * m) / denom if denom > 0 else 1.0
    d0 = float(np.clip(d0, 1e-2, 1e6))
    s02 = m * df_res / (df_res + d0)          # 先验方差(收缩均值)
    s2_post = (df_res * sigma2 + d0 * s02) / (df_res + d0)
    var_grp = np.diag(XtX_inv)[1] * s2_post  # grp 系数方差
    se_mod = np.sqrt(var_grp)
    se_mod = np.clip(se_mod, 1e-4, None)      # 防 0 权重
    fc = beta[1, :]
    modt = fc / se_mod
    df_tot = df_res + d0
    pval = 2 * stats.t.sf(np.abs(modt), df_tot)
    out = pd.DataFrame({
        "log2FC": fc,
        "modt": modt,
        "P.Value": pval,
        "SE": se_mod,
    }, index=expr.index)
    return out

def bh(p):
    p = np.asarray(p, dtype=float)
    p = np.where(np.isnan(p), 1.0, p)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    return out

def deg_table(limma_res, fc_thr=0.5, pthr=0.05):
    res = limma_res.copy()
    res["adj.P.Val"] = bh(res["P.Value"].values)
    res["direction"] = np.where(res["log2FC"] > 0, "up", "down")
    sig = res[(res["log2FC"].abs() > fc_thr) & (res["adj.P.Val"] < pthr)]
    return res, sig

# ----------------------------------------------------------------------
# 2. 各数据集分组定义
# ----------------------------------------------------------------------
def group_gse98895(meta):
    titles = meta["!Sample_title"]
    grp = []
    for t in titles:
        t = t.strip('"')
        suf = t.split("-")[-1].strip().upper()
        grp.append(1 if suf in ("MS", "SM") else 0)   # MS+SM=病例, CTRL=对照
    return np.array(grp), grp.count(1), len(grp) - grp.count(1)

def group_gse21510(meta):
    chars = meta["!Sample_characteristics_ch1"]
    grp = []
    for c in chars:
        c = c.strip('"').lower()
        if c.startswith("tissue: cancer"):
            grp.append(1)
        elif c.startswith("tissue: normal"):
            grp.append(0)
        else:
            grp.append(-1)   # 丢弃
    return np.array(grp)

def group_gse32323(meta):
    chars = meta["!Sample_characteristics_ch1"]
    grp = []
    for c in chars:
        c = c.strip('"').lower()
        if c == "tissue: cancer, homogenized":
            grp.append(1)
        elif c == "tissue: normal, homogenized":
            grp.append(0)
        else:
            grp.append(-1)   # 细胞系等丢弃
    return np.array(grp)

def group_gse39582(meta):
    src = meta["!Sample_source_name_ch1"]
    grp = []
    for s in src:
        s = s.strip('"').lower()
        if "adenocarcinoma" in s:
            grp.append(1)
        elif "non tumoral" in s:
            grp.append(0)
        else:
            grp.append(-1)
    return np.array(grp)

# ----------------------------------------------------------------------
# 3. 主流程
# ----------------------------------------------------------------------
def main():
    sym_map = load_annot(ANNOT_GPL6947)
    print(f"[annot] GSE98895 真实平台 GPL6947(Illumina) 探针->基因映射数 = {len(sym_map)}")

    # ============ Stage 1: GSE98895 ============
    print("\n========== Stage 1: GSE98895 (MetS vs 对照) ==========")
    p1 = os.path.join(BASE, "GSE98895", "GSE98895_series_matrix.txt.gz")
    meta = parse_meta(p1)
    grp, ncase, nctrl = group_gse98895(meta)
    print(f"  分组: MetS(MS+SM)={ncase}, 对照(CTRL)={nctrl}")
    expr = read_expr(p1)
    expr = probes_to_genes(expr, sym_map)
    expr = maybe_log2(expr)
    print(f"  基因级表达矩阵: {expr.shape}")
    lm = limma_modt(expr, grp)
    res, sig = deg_table(lm)
    res.to_csv(os.path.join(BASE, "GSE98895", "DEG_MetS_vs_Control.csv"))
    sig.to_csv(os.path.join(BASE, "GSE98895", "DEG_MetS_vs_Control.sig.csv"))
    print(f"  DEG (|log2FC|>0.5 & adj<0.05): {len(sig)} 个")
    print("  Top10:")
    print(sig.reindex(sig["P.Value"].abs().sort_values().index).head(10)[["log2FC","P.Value","adj.P.Val","direction"]])

    # ============ Stage 2: MetS 血液基因级 IVW meta (GSE98895 Illumina + GSE200744 RNA-seq) ============
    print("\n========== Stage 2: MetS 血液基因级 IVW meta (GSE98895 + GSE200744) ==========")
    # Study A: GSE98895 limma (gene-level), 已算于 Stage1 -> lm
    fcA = lm["log2FC"].copy()
    seA = lm["SE"].copy()
    # Study B: GSE200744 DESeq2 全基因结果 (log2FC + lfcSE)
    g200_full = pd.read_csv(os.path.join(BASE, "GSE200744", "DEG_MetS_vs_Normal.csv"))
    g200_full = g200_full.dropna(subset=["Gene Symbol", "log2FC", "lfcSE"])
    g200_full = g200_full[g200_full["Gene Symbol"].astype(str) != ""]
    g200_full = g200_full.groupby("Gene Symbol").first()
    fcB = g200_full["log2FC"]; seB = g200_full["lfcSE"]
    common_m = sorted(set(fcA.index) & set(fcB.index))
    print(f"  MetS 两研究共有基因 = {len(common_m)}")
    Fc = pd.DataFrame({"GSE98895": fcA.reindex(common_m), "GSE200744": fcB.reindex(common_m)})
    Se = pd.DataFrame({"GSE98895": seA.reindex(common_m), "GSE200744": seB.reindex(common_m)})
    w = 1.0 / (Se ** 2)
    mfc = (Fc * w).sum(axis=1) / w.sum(axis=1)
    mse = np.sqrt(1.0 / w.sum(axis=1))
    mz = mfc / mse
    mp = 2 * stats.norm.sf(mz.abs())
    mres = pd.DataFrame({"log2FC": mfc, "SE": mse, "Z": mz, "P.Value": mp})
    mres["adj.P.Val"] = bh(mp)
    mres["direction"] = np.where(mres["log2FC"] > 0, "up", "down")
    mres = mres.sort_values("P.Value")
    mres.to_csv(os.path.join(RES, "MetS_meta_DGE.csv"))
    mets_sig = mres[(mres["log2FC"].abs() > 0.5) & (mres["adj.P.Val"] < 0.05)]
    mets_sig.to_csv(os.path.join(RES, "MetS_meta_DGE.sig.csv"))
    print(f"  MetS meta DEG (|fc|>0.5 & adj<0.05): {len(mets_sig)} 个")
    print("  Top10:")
    print(mets_sig.head(10)[["log2FC", "P.Value", "adj.P.Val", "direction"]])
    mets_fc = mres["log2FC"].copy()   # 用于 Stage4 方向判断

    # ============ Stage 3: CRC 训练集 GPL570 x3 meta ============
    print("\n========== Stage 3: CRC 训练集 GPL570 x3 meta ==========")
    crc_studies = {
        "GSE21510": (group_gse21510, 148),
        "GSE32323": (group_gse32323, 44),
        "GSE39582": (group_gse39582, 585),
    }
    sym_map_crc = load_annot(ANNOT_GPL570)
    print(f"[annot] CRC 训练集 GPL570 探针->基因映射数 = {len(sym_map_crc)}")
    study_fc, study_se = {}, {}
    for gse, (gfunc, expected) in crc_studies.items():
        print(f"  --- {gse} ---")
        pp = os.path.join(BASE, gse, f"{gse}_series_matrix.txt.gz")
        meta = parse_meta(pp)
        grp = gfunc(meta)
        keep = grp >= 0
        ncase = int((grp[keep] == 1).sum()); nctrl = int((grp[keep] == 0).sum())
        print(f"    分组: tumor={ncase}, normal={nctrl} (丢弃={(~keep).sum()})")
        expr = read_expr(pp)
        expr = probes_to_genes(expr, sym_map_crc)
        expr = maybe_log2(expr)
        expr = expr.loc[:, keep]                 # 仅保留肿瘤/正常
        grp_k = grp[keep]
        lm = limma_modt(expr, grp_k)
        # 仅保留两组的基因(都应有估计); 存 fc & se
        fc = lm["log2FC"]; se = lm["SE"]
        study_fc[gse] = fc
        study_se[gse] = se
        lm_out = pd.DataFrame({"log2FC": fc, "SE": se, "P.Value": lm["P.Value"]})
        lm_out.to_csv(os.path.join(BASE, gse, f"DEG_tumor_vs_normal.csv"))
        print(f"    基因数={expr.shape[0]}")

    # IVW meta: 取三研究共有基因
    common_genes = set(study_fc["GSE21510"].index)
    for gse in ("GSE32323", "GSE39582"):
        common_genes &= set(study_fc[gse].index)
    common_genes = sorted(common_genes)
    print(f"  三研究共有基因数 = {len(common_genes)}")
    Fc = pd.DataFrame({g: study_fc[g].reindex(common_genes) for g in crc_studies})
    Se = pd.DataFrame({g: study_se[g].reindex(common_genes) for g in crc_studies})
    w = 1.0 / (Se ** 2)
    meta_fc = (Fc * w).sum(axis=1) / w.sum(axis=1)
    meta_se = np.sqrt(1.0 / w.sum(axis=1))
    meta_z = meta_fc / meta_se
    meta_p = 2 * stats.norm.sf(meta_z.abs())
    meta_res = pd.DataFrame({
        "log2FC": meta_fc, "SE": meta_se, "Z": meta_z, "P.Value": meta_p,
    })
    meta_res["adj.P.Val"] = bh(meta_p)
    meta_res["direction"] = np.where(meta_res["log2FC"] > 0, "up", "down")
    meta_res = meta_res.sort_values("P.Value")
    meta_res.to_csv(os.path.join(RES, "CRC_meta_DGE.csv"))
    crc_sig = meta_res[(meta_res["log2FC"].abs() > 0.5) & (meta_res["adj.P.Val"] < 0.05)]
    crc_sig.to_csv(os.path.join(RES, "CRC_meta_DGE.sig.csv"))
    print(f"  CRC meta DEG (|fc|>0.5 & adj<0.05): {len(crc_sig)} 个")
    print("  Top10:")
    print(crc_sig.head(10)[["log2FC","P.Value","adj.P.Val","direction"]])

    # ============ Stage 4: co-DEGs (MetS meta sig ∩ CRC meta sig, 同方向) ============
    print("\n========== Stage 4: 跨病 co-DEGs ==========")
    crc_sig_fc = crc_sig["log2FC"].copy()
    co = []
    for g in mets_sig.index:
        if g in crc_sig_fc.index and np.sign(mets_fc[g]) == np.sign(crc_sig_fc[g]):
            co.append(g)
    print(f"  MetS meta sig={len(mets_sig)}  CRC meta sig={len(crc_sig)}  co-DEGs(同方向)={len(co)}")
    print("  ", co)
    co_df = pd.DataFrame({"Gene": co})
    co_df["MetS_meta_log2FC"] = mets_fc.reindex(co).values
    co_df["CRC_meta_log2FC"] = crc_sig_fc.reindex(co).values
    co_df.to_csv(os.path.join(RES, "co_DEGs.csv"), index=False)

    print("\n========== 完成 ==========")
    print(f"输出目录: {RES}")
    print("  GSE98895/DEG_MetS_vs_Control.csv(.sig.csv)")
    print("  results/MetS_meta_DGE.csv(.sig.csv)   <- MetS 血液共识(IVW)")
    print("  results/CRC_meta_DGE.csv(.sig.csv)")
    print("  results/co_DEGs.csv")

if __name__ == "__main__":
    main()

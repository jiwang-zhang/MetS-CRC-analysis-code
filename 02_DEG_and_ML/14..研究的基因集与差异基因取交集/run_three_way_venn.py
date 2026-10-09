"""
三集合交集 + 漂亮 Venn 图
集合: CRC差异基因(4319) ∩ MetS差异基因(830) ∩ 研究基因集 TNF/NF-kB freq2(1803)
输出: 各区域基因列表CSV + 高分辨率Venn图(PNG/PDF)
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

BASE = r"D:\daixie-colorectal\第二部分DEG"
F11  = os.path.join(BASE, "11..研究的基因集与差异基因取交集")
OUT  = F11
os.makedirs(OUT, exist_ok=True)

def load_genes(path, col=0):
    """读取基因列表（第一列或指定列），转大写集合"""
    genes = set()
    with open(path, encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, None)  # 跳过表头
        for row in reader:
            if not row:
                continue
            g = row[col].strip().upper()
            if g:
                genes.add(g)
    return genes

# ---------- 1. 载入三集合 ----------
crc  = load_genes(os.path.join(BASE, "09.两个疾病分别DEG", "DE_significant_CRC.csv"))
mets = load_genes(os.path.join(BASE, "09.两个疾病分别DEG", "DE_significant_MetS.csv"))
rese = load_genes(os.path.join(F11, "TNF_NFKB_freq2_with_counts.csv"))

print(f"CRC 差异基因   : {len(crc)}")
print(f"MetS 差异基因  : {len(mets)}")
print(f"研究基因集(TNF/NFKB freq2): {len(rese)}")

# ---------- 2. 7 个区域 ----------
only_crc   = crc - mets - rese
only_mets  = mets - crc - rese
only_rese  = rese - crc - mets
crc_mets   = (crc & mets) - rese          # 共差异基因(不含研究集)
crc_rese   = (crc & rese) - mets
mets_rese  = (mets & rese) - crc
triple     = crc & mets & rese            # ★ 三交集(核心结果)

regions = {
    "only_CRC": only_crc,
    "only_MetS": only_mets,
    "only_Research": only_rese,
    "CRC_and_MetS_notResearch": crc_mets,
    "CRC_and_Research_notMetS": crc_rese,
    "MetS_and_Research_notCRC": mets_rese,
    "CRC_and_MetS_and_Research": triple,
}

print("\n===== 7 个区域基因数 =====")
for k, v in regions.items():
    print(f"  {k:32s}: {len(v)}")

# ---------- 3. 导出各区域基因CSV ----------
def dump(name, genes):
    p = os.path.join(OUT, f"region_{name}.csv")
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Gene"])
        for g in sorted(genes):
            w.writerow([g])
    return p

for k, v in regions.items():
    dump(k, v)

# 三交集附带来源标注(研究集里的founder计数)
rese_count = {}
with open(os.path.join(F11, "TNF_NFKB_freq2_with_counts.csv"), encoding="utf-8") as f:
    r = csv.DictReader(f)
    for row in r:
        rese_count[row["Gene"].strip().upper()] = row.get("FounderSetCount", "")

with open(os.path.join(OUT, "region_triple_intersection_annotated.csv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Gene", "TNF_NFKB_FounderSetCount"])
    for g in sorted(triple):
        w.writerow([g, rese_count.get(g, "")])

print(f"\n★ 三交集(CRC∩MetS∩TNF/NF-kB)基因数 = {len(triple)}")
print("  基因:", ", ".join(sorted(triple)))

# ---------- 4. 漂亮 Venn 图 ----------
plt.rcParams["font.family"] = "DejaVu Sans"
fig, ax = plt.subplots(figsize=(9, 7.5), dpi=200)

# 颜色: CRC=青绿, MetS=橙, Research=紫
v = venn3(
    [crc, mets, rese],
    set_labels=("CRC DEGs\n(4,319)", "MetS DEGs\n(830)", "TNF/NF-κB\nsignature (1,803)"),
    set_colors=("#2ca19f", "#e8893b", "#7b5ea7"),
    alpha=0.55,
    ax=ax,
)
# 数字标签美化
if v:
    for idv in v.subset_labels:
        if idv:
            idv.set_fontsize(12)
            idv.set_fontweight("bold")
            idv.set_color("#222222")
    for idv in v.set_labels:
        if idv:
            idv.set_fontsize(13)
            idv.set_fontweight("bold")
    # 加描边让圆更清爽
    c = venn3_circles([crc, mets, rese], linestyle="solid", linewidth=1.4, color="#444444", ax=ax)

ax.set_title(
    "Three-way intersection of CRC DEGs, MetS DEGs and TNF/NF-κB signature",
    fontsize=14, fontweight="bold", pad=14, color="#1a1a1a",
)
# 副标题标注三交集核心数
ax.text(0.5, -0.62,
        f"Core intersection (CRC ∩ MetS ∩ TNF/NF-κB) = {len(triple)} genes",
        ha="center", fontsize=11, fontstyle="italic", color="#7b5ea7",
        transform=ax.transAxes)

plt.tight_layout()
png = os.path.join(OUT, "three_way_venn.png")
pdf = os.path.join(OUT, "three_way_venn.pdf")
plt.savefig(png, dpi=300, bbox_inches="tight")
plt.savefig(pdf, bbox_inches="tight")
plt.close()
print(f"\nVenn 图已保存:\n  {png}\n  {pdf}")

# 汇总表
with open(os.path.join(OUT, "intersection_summary.csv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Region", "GeneCount"])
    for k, v in regions.items():
        w.writerow([k, len(v)])
    w.writerow(["TOTAL_unique", len(crc | mets | rese)])
print("\n全部完成。")

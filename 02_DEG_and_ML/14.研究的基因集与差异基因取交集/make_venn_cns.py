"""
CNS / 顶刊出版级 三集合 Venn 图
集合: CRC DEGs(4319) ∩ MetS DEGs(830) ∩ TNF/NF-κB signature(1803)
风格: Nature 高级配色 + 深色描边 + 白色数字(深色描边) + 柔和投影 + 600 dpi
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import matplotlib.patheffects as pe

BASE = r"D:\daixie-colorectal\第二部分DEG"
F11  = os.path.join(BASE, "11..研究的基因集与差异基因取交集")
OUT  = F11

def load_genes(path, col=0):
    s = set()
    with open(path, encoding="utf-8") as f:
        r = csv.reader(f)
        next(r, None)
        for row in r:
            if row:
                g = row[col].strip().upper()
                if g:
                    s.add(g)
    return s

crc  = load_genes(os.path.join(BASE, "09.两个疾病分别DEG", "DE_significant_CRC.csv"))
mets = load_genes(os.path.join(BASE, "09.两个疾病分别DEG", "DE_significant_MetS.csv"))
rese = load_genes(os.path.join(F11, "TNF_NFKB_freq2_with_counts.csv"))

only_crc  = crc - mets - rese
only_mets = mets - crc - rese
only_rese = rese - crc - mets
crc_mets  = (crc & mets) - rese
crc_rese  = (crc & rese) - mets
mets_rese = (mets & rese) - crc
triple    = crc & mets & rese

# ---- Nature 风格高级配色 ----
C_CRC, C_METS, C_RESE = "#1F7A8C", "#E3812B", "#B0386B"   # fill
E_CRC, E_METS, E_RESE = "#0C4651", "#8F4A12", "#6E2142"   # 深描边同系
EDGES = [E_CRC, E_METS, E_RESE]

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.weight'] = 'normal'

fig, ax = plt.subplots(figsize=(9.6, 8.6), dpi=600)
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
ax.set_axis_off()

v = venn3(
    [crc, mets, rese],
    set_labels=("CRC DEGs  (4,319)", "MetS DEGs  (830)", "TNF/NF-\u03baB signature  (1,803)"),
    set_colors=(C_CRC, C_METS, C_RESE),
    alpha=0.62,
    ax=ax,
)

if v:
    # 区域数字: 白字 + 深色描边(任意背景都清晰)
    for lbl in v.subset_labels:
        if lbl:
            lbl.set_fontsize(13)
            lbl.set_fontweight('bold')
            lbl.set_color('white')
            lbl.set_path_effects([pe.withStroke(linewidth=3.0, foreground='#262626')])
    # 集合标签: 外侧粗体深色
    for lbl in v.set_labels:
        if lbl:
            lbl.set_fontsize(15)
            lbl.set_fontweight('bold')
            lbl.set_color('#1a1a1a')
    # 椭圆: 深描边 + 柔和投影
    for patch, edge in zip(v.patches, EDGES):
        patch.set_edgecolor(edge)
        patch.set_linewidth(2.0)
        patch.set_path_effects([
            pe.SimpleLineShadow(offset=(2.2, -2.2), alpha=0.16, shadow_color='#555555'),
            pe.Normal(),
        ])

# 标题与脚注(放在 figure 级, 避免被裁切)
fig.text(0.5, 0.955,
         "Three-way intersection of CRC and MetS DEGs with the TNF/NF-\u03baB signature",
         ha='center', va='center', fontsize=14, fontweight='bold', color='#1a1a1a')
fig.text(0.5, 0.045,
         f"Core intersection (CRC \u2229 MetS \u2229 TNF/NF-\u03baB) = {len(triple)} genes",
         ha='center', va='center', fontsize=12, fontstyle='italic', color=C_RESE)

plt.subplots_adjust(top=0.90, bottom=0.10, left=0.06, right=0.94)
png = os.path.join(OUT, "venn_CNS.png")
pdf = os.path.join(OUT, "venn_CNS.pdf")
plt.savefig(png, dpi=600, bbox_inches="tight", facecolor='white')
plt.savefig(pdf, bbox_inches="tight", facecolor='white')
plt.close()
print(f"Venn(CNS) saved:\n  {png}\n  {pdf}")
print(f"Triple intersection = {len(triple)} genes")

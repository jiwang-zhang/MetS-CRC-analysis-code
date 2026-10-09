import os, gzip
import numpy as np, pandas as pd
BASE = r"D:/daixie-colorectal/19.databases"
p = os.path.join(BASE, "GSE98895", "GSE98895_series_matrix.txt.gz")

# 1) read_expr 行为
expr = pd.read_csv(p, sep="\t", index_col=0, comment="!", compression="gzip", low_memory=False)
print("read_expr shape:", expr.shape)
print("dtypes sample:\n", expr.dtypes.head(3))
print("index head:", list(expr.index[:3]))
print("values head:\n", expr.iloc[:3, :3])
print("first cell raw repr:", repr(expr.iloc[0,0]), type(expr.iloc[0,0]))

# 2) annot 映射检查
sym_map = {}
with gzip.open(os.path.join(BASE,"GSE98895","GPL570.annot.gz"),"rt",errors="ignore") as f:
    header=None
    for line in f:
        if line.startswith("!platform_table_begin"): continue
        if line.startswith("#") or line.startswith("!") or line.startswith("^"): continue
        if header is None:
            header=line.rstrip("\n").split("\t"); i_id=header.index("ID"); i_sym=header.index("Gene symbol"); continue
        parts=line.rstrip("\n").split("\t")
        if len(parts)<=max(i_id,i_sym): continue
        pid=parts[i_id]; sym=parts[i_sym].split("///")[0].strip()
        if sym and sym not in ("","---","NA","NaN"): sym_map[pid]=sym
print("annot map size:", len(sym_map))
print("expr index sample -> sym_map:", [(i, sym_map.get(str(i),None)) for i in expr.index[:5]])

# 3) 模拟 probes_to_genes
idx = [sym_map.get(str(x), None) for x in expr.index]
print("mapped None count:", sum(1 for x in idx if x is None), "/", len(idx))
expr2 = expr.copy()
expr2.index = idx
expr2 = expr2[~expr2.index.isna()]
print("after index filter shape:", expr2.shape)
expr2.index = expr2.index.astype(str)
expr2 = expr2.apply(pd.to_numeric, errors="coerce")
print("after to_numeric dtypes head:\n", expr2.dtypes.head(3))
print("any finite?:", np.isfinite(expr2.to_numpy()).any())
print("sample values:\n", expr2.iloc[:3,:3])

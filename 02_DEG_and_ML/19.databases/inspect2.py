import gzip, os
BASE = r"D:/daixie-colorectal/19.databases"

# 1) GSE98895 all titles
p = os.path.join(BASE, "GSE98895", "GSE98895_series_matrix.txt.gz")
with gzip.open(p, "rt", errors="ignore") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_title"):
            titles = line.split("\t")[1:]
suffs = {}
for t in titles:
    t = t.strip('"')
    suf = t.split(" - ")[-1] if " - " in t else t
    suffs[suf] = suffs.get(suf, 0) + 1
print("GSE98895 title suffixes:", suffs)
print("All 40 titles:")
for t in titles:
    print("  ", t.strip('"'))

# 2) GPL570.annot column header
print("\nGPL570.annot.gz structure:")
with gzip.open(os.path.join(BASE, "GSE98895", "GPL570.annot.gz"), "rt", errors="ignore") as f:
    lines = f.readlines()
for i, l in enumerate(lines[:60]):
    if "Gene Symbol" in l or l.startswith("ID\t") or "IDENTIFIER" in l.upper():
        print(f"  [line {i}] {l.rstrip(chr(10))[:200]}")
        # print next 2 rows
        for j in range(i+1, min(i+3, len(lines))):
            print(f"  [line {j}] {lines[j].rstrip(chr(10))[:200]}")
        break

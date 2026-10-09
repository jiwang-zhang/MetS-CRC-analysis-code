import gzip, os
BASE = r"D:/daixie-colorectal/19.databases"
with gzip.open(os.path.join(BASE, "GSE98895", "GPL570.annot.gz"), "rt", errors="ignore") as f:
    header = None
    for line in f:
        if line.startswith("#") or line.startswith("!") or line.startswith("^") or line.startswith(">"):
            continue
        # first data-like line with tabs and 'Gene Symbol'
        if "\t" in line and ("Gene Symbol" in line or line.startswith("ID")):
            header = line.rstrip("\n").split("\t")
            break
print("NCOLS:", len(header))
for i, c in enumerate(header):
    if c in ("ID", "Gene Symbol", "ENTREZ_GENE_ID", "Gene Title"):
        print(f"  col {i}: {c}")
# print 3 sample rows' ID + Gene Symbol
idx_id = header.index("ID")
idx_sym = header.index("Gene Symbol")
with gzip.open(os.path.join(BASE, "GSE98895", "GPL570.annot.gz"), "rt", errors="ignore") as f:
    cnt = 0
    for line in f:
        if line.startswith("#") or line.startswith("!") or line.startswith("^") or line.startswith(">"):
            continue
        if "\t" not in line:
            continue
        parts = line.rstrip("\n").split("\t")
        if parts[0] == "ID":
            continue
        print(f"  {parts[idx_id]} -> {parts[idx_sym]}")
        cnt += 1
        if cnt >= 5:
            break

# -*- coding: utf-8 -*-
"""Open Targets Platform v4 GraphQL: 取疾病-靶点(基因)关联基因集，与 co-DEG 取交集。"""
import urllib.request, json, os

BASE = r"D:\daixie-colorectal\第二部分DEG"
CO_DEG = os.path.join(BASE, r"09.两个疾病分别DEG", "co_DEGs_CRC_MetS.csv")
OUT = r"D:\daixie-colorectal\第二部分DEG\10.获取研究的基因集\hot_databases"
URL = "https://api.platform.opentargets.org/api/v4/graphql"

def gql(query):
    # Open Targets 要求匿名查询外层包 {}
    body = "{" + query + "}"
    req = urllib.request.Request(URL, data=json.dumps({"query": body}).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())

# 1) 搜疾病 ID（不带 entity 过滤，python 内筛选）
diseases = {
    "ColorectalCancer": "colorectal cancer",
    "Type2Diabetes": "type 2 diabetes",
    "Obesity": "obesity",
    "MetabolicSyndrome": "metabolic syndrome",
}
ids = {}
for key, q in diseases.items():
    res = gql(f'search(queryString:"{q}"){{ hits{{ id name entity }} }}')
    hits = res.get("data", {}).get("search", {}).get("hits", [])
    dis_hits = [h for h in hits if h.get("entity") == "disease"]
    if dis_hits:
        ids[key] = (dis_hits[0]["id"], dis_hits[0]["name"])
        print(f"[search] {key}: {dis_hits[0]['id']} | {dis_hits[0]['name']}")
    else:
        print(f"[search] {key}: NO disease hit (got {[h.get('entity') for h in hits[:3]]})")

# 2) 取每个疾病的关联靶点基因（Open Targets 单页上限 3000，需翻页）
disease_genes = {}
for key, (did, dname) in ids.items():
    syms = []
    idx = 0
    while True:
        q = f'disease(efoId:"{did}"){{ associatedTargets(page:{{index:{idx},size:3000}}){{ count rows{{ target{{ approvedSymbol }} }} }} }}'
        res = gql(q)
        dis = res.get("data", {}).get("disease")
        if not dis or not dis.get("associatedTargets"):
            break
        rows = dis["associatedTargets"]["rows"]
        if not rows:
            break
        syms += [r["target"]["approvedSymbol"].upper() for r in rows]
        total = dis["associatedTargets"]["count"]
        if len(syms) >= total or len(rows) < 3000:
            break
        idx += 1
    if syms:
        disease_genes[key] = syms
        print(f"[targets] {key} ({dname}): {len(syms)} genes (共 {total} 关联靶点)")
        with open(os.path.join(OUT, f"ot_{key}.txt"), "w", encoding="utf-8") as f:
            f.write("\n".join(sorted(set(syms))))
    else:
        print(f"[targets] {key}: null/empty")

# 3) 与 200 co-DEG 取交集
co = []
with open(CO_DEG, encoding="utf-8") as f:
    next(f)
    for line in f:
        g = line.split(",")[0].strip().strip('"')
        if g:
            co.append(g.upper())
co_set = set(co)
print(f"\n[co-DEG] {len(co_set)} genes")

with open(os.path.join(OUT, "opentargets_intersections.csv"), "w", encoding="utf-8") as f:
    f.write("disease,ot_gene_count,n_hits,hit_genes\n")
    for key, syms in disease_genes.items():
        sset = set(syms)
        hit = sorted(co_set & sset)
        print(f"\n=== {key}: {len(hit)} co-DEG 命中 ===")
        print(", ".join(hit) if hit else "(无)")
        f.write(f"{key},{len(sset)},{len(hit)},\";\".join(hit)\n")

# -*- coding: utf-8 -*-
"""检索 TNF-a / NF-kB 轴引出 HALLMARK_TNFA_SIGNALING_VIA_NFKB 的支撑文献。
E-utilities (esearch + esummary)，标准库，无需第三方包。
"""
import urllib.request, urllib.parse, json, time, sys

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

QUERIES = {
    "CRC_TNF_NFkB": '(colorectal cancer[tiab] OR colorectal neoplasms[mh]) AND (NF-kappaB[tiab] OR NF-kB[tiab] OR TNF-alpha[tiab] OR "tumor necrosis factor"[tiab]) AND (inflammation[tiab] OR tumorigenesis[tiab] OR review[pt]) AND 2015:2026[dp]',
    "MetS_TNF_NFkB": '(metabolic syndrome[tiab] OR obesity[mh] OR insulin resistance[tiab]) AND (NF-kappaB[tiab] OR NF-kB[tiab] OR TNF-alpha[tiab] OR "tumor necrosis factor"[tiab]) AND (chronic inflammation[tiab] OR "low-grade inflammation"[tiab]) AND 2015:2026[dp]',
    "MetS_CRC_bridge": '(metabolic syndrome[tiab] OR obesity[tiab] OR insulin resistance[tiab]) AND (colorectal cancer[tiab] OR colorectal neoplasms[mh]) AND (inflammation[tiab] OR NF-kappaB[tiab] OR TNF-alpha[tiab]) AND 2015:2026[dp]',
}

def esearch(term, retmax=8):
    params = urllib.parse.urlencode({
        "db": "pubmed", "term": term, "retmax": retmax,
        "retmode": "json", "sort": "relevance"})
    url = BASE + "esearch.fcgi?" + params
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.load(r)["esearchresult"]["idlist"]

def esummary(pmids):
    params = urllib.parse.urlencode({
        "db": "pubmed", "id": ",".join(pmids), "retmode": "json"})
    url = BASE + "esummary.fcgi?" + params
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.load(r)["result"]

out = []
for tag, term in QUERIES.items():
    out.append("\n" + "=" * 70)
    out.append(f"[{tag}]")
    out.append(term)
    out.append("-" * 70)
    try:
        ids = esearch(term)
        time.sleep(0.4)
        if not ids:
            out.append("  (no hits)")
            continue
        res = esummary(ids)
        time.sleep(0.4)
        for pid in ids:
            it = res.get(pid, {})
            title = it.get("title", "")
            journal = it.get("source", "")
            year = (it.get("pubdate", "") or "")[:4]
            out.append(f"  PMID:{pid} | {journal} {year}")
            out.append(f"    {title}")
    except Exception as e:
        out.append(f"  ERROR: {e}")

text = "\n".join(out)
print(text)
with open("lit_tnfa_nfkb_hits.txt", "w", encoding="utf-8") as f:
    f.write(text)

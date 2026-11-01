# ======================================================================
# asn_build_records.py
#
# Purpose  : Build one citation record per ASN 2016-18 application from the frozen OpenAlex pull:
#            the CV works published in the ASN window, with citation counts reconstructed to the
#            application deadline. Three versions of the deadline-year citations are kept:
#            "prorata" (primary; the deadline year's citations are prorated by the elapsed fraction
#            of the year), "lower" (none of the deadline year's citations) and "upper" (all of them).
# Reads    : data/asn/oa_raw/selection.tsv, data/asn/oa_raw/pages.jsonl, data/asn/oa_raw/single.jsonl,
#            data/asn/derived/asn1618_bibl_linked_noNames.tsv (DOI lists)
# Produces : data/asn/derived/asn_records.jsonl  one line per application
#            data/asn/derived/asn_records_summary.tsv  per-application summary statistics
# Run      : from the repository root -- python R/asn_build_records.py
# Notes    : Window: works published from 1 January of year T - 15 (full professor, level 1) or
#            T - 10 (associate professor, level 2) up to the deadline, T = deadline year; this is the
#            ASN 2016 window for the citation and h-index indicators. Works are de-duplicated by
#            OpenAlex id. The deadlines are the application deadlines of the five 2016-18 terms
#            (D.D. 29 July 2016 n. 1532, art. 2), at which the official indicators were computed.
# ======================================================================
import ast, json, os, datetime as dt
import pandas as pd

R = os.path.join("data", "asn")
DEAD = {"1": dt.date(2016, 12, 2), "2": dt.date(2017, 4, 3), "3": dt.date(2017, 8, 4),
        "4": dt.date(2017, 12, 5), "5": dt.date(2018, 4, 6)}   # application deadlines, D.D. 1532/2016 art. 2 c. 1

W = {}
for line in open(os.path.join(R, "oa_raw", "pages.jsonl"), encoding="utf-8"):
    for w in json.loads(line)["results"] or []:
        if w.get("doi"):
            d = w["doi"].lower().replace("https://doi.org/", "")
            if d not in W or (w.get("cited_by_count") or 0) > (W[d].get("cited_by_count") or 0):
                W[d] = w
for line in open(os.path.join(R, "oa_raw", "single.jsonl"), encoding="utf-8"):
    j = json.loads(line)
    if j["result"] and "id" in j["result"]:
        W.setdefault(j["doi"], j["result"])

sel = pd.read_csv(os.path.join(R, "oa_raw", "selection.tsv"), sep="\t", dtype=str)
L = pd.read_csv(os.path.join(R, "derived", "asn1618_bibl_linked_noNames.tsv"), sep="\t", dtype=str)[["id", "dois"]]
sel = sel.merge(L, on="id", how="left")

def asof(w, d, mode):
    cby = {c["year"]: c["cited_by_count"] for c in (w.get("counts_by_year") or [])}
    later = sum(v for y, v in cby.items() if y > d.year)
    cur = cby.get(d.year, 0)
    base = (w.get("cited_by_count") or 0) - later - cur
    if mode == "lower":
        x = base
    elif mode == "upper":
        x = base + cur
    else:
        f = (d - dt.date(d.year, 1, 1)).days / 365.25
        x = base + f * cur
    return max(x, 0.0)

def hindex(c):
    c = sorted(c, reverse=True)
    return sum(1 for i, x in enumerate(c) if x >= i + 1)

out = open(os.path.join(R, "derived", "asn_records.jsonl"), "w", encoding="utf-8")
rows = []
for _, p in sel.iterrows():
    d = DEAD[p.quadrimestre]
    start = d.year - (15 if p.fascia == "1" else 10)
    ws, seen = [], set()
    for x in ast.literal_eval(p.dois):
        w = W.get(x.strip().lower())
        if w and w["id"] not in seen:
            seen.add(w["id"]); ws.append(w)
    win = [w for w in ws if w.get("publication_year") is not None and start <= w["publication_year"] <= d.year
           and (w.get("publication_date") or "9999") <= d.isoformat()]
    rec = {"id": p.id, "term": p.quadrimestre, "level": p.fascia, "sector": p.settore, "deadline": d.isoformat(),
           "I1": int(float(p.I1)), "I2": int(float(p.I2)), "I3": int(float(p.I3)),
           "S1": int(float(p.S1_art)), "S2": int(float(p.S2_cit)), "S3": int(float(p.S3_h)),
           "n_doi": int(p.n_doi), "n_matched": len(ws), "n_window": len(win)}
    for mode in ("prorata", "lower", "upper"):
        rec["x_" + mode] = [round(asof(w, d, mode), 3) for w in win]
    out.write(json.dumps(rec) + "\n")
    xs = sorted(rec["x_prorata"])
    n = len(xs)
    xi = [int(v) for v in rec["x_prorata"]]
    rows.append({k: rec[k] for k in ("id", "term", "level", "sector", "I3", "S3", "n_doi", "n_matched", "n_window")} |
                {"H_prorata": hindex(xi), "H_lower": hindex([int(v) for v in rec["x_lower"]]),
                 "H_upper": hindex([int(v) for v in rec["x_upper"]]),
                 "p0": sum(1 for v in xi if v == 0) / n if n else None,
                 "q50": pd.Series(xi).quantile(0.5) if n else None, "q90": pd.Series(xi).quantile(0.9) if n else None})
out.close()
pd.DataFrame(rows).to_csv(os.path.join(R, "derived", "asn_records_summary.tsv"), sep="\t", index=False)
print("applications", len(rows))

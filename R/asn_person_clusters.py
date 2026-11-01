# ======================================================================
# asn_person_clusters.py
#
# Purpose  : Group ASN applications that are likely to belong to the same person (applications in
#            several sectors, levels or terms), using overlap of the CV DOI lists only (no names).
#            Two applications are linked when the Jaccard overlap of their DOI sets is at least 0.5;
#            clusters are the connected components. Used to exclude a candidate's own other
#            applications from the donors of pooled estimators and to resample candidates, not
#            applications, in the bootstrap.
# Reads    : data/asn/derived/asn1618_bibl_linked_noNames.tsv, data/asn/oa_raw/selection.tsv
# Produces : data/asn/derived/asn_person_clusters.tsv  (id, cluster, cluster_size)
# Run      : from the repository root -- python R/asn_person_clusters.py
# ======================================================================
import ast, os, collections
import pandas as pd

R = os.path.join("data", "asn")
L = pd.read_csv(os.path.join(R, "derived", "asn1618_bibl_linked_noNames.tsv"), sep="\t", dtype=str)[["id", "dois"]]
sel = pd.read_csv(os.path.join(R, "oa_raw", "selection.tsv"), sep="\t", dtype=str)[["id"]]
L = L[L.id.isin(set(L.id))]  # all linked applications take part, so clusters span the selection and its duplicates
D = {r.id: {d.strip().lower() for d in ast.literal_eval(r.dois) if d.strip()} for r in L.itertuples()}

inv = collections.defaultdict(list)
for a, ds in D.items():
    for d in ds:
        inv[d].append(a)
parent = {a: a for a in D}
def find(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]; a = parent[a]
    return a
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra
checked = set()
for d, apps in inv.items():
    if len(apps) < 2 or len(apps) > 50:
        continue
    for i in range(len(apps)):
        for j in range(i + 1, len(apps)):
            a, b = apps[i], apps[j]
            key = (a, b) if a < b else (b, a)
            if key in checked:
                continue
            checked.add(key)
            A, B = D[a], D[b]
            if A and B and len(A & B) / len(A | B) >= 0.5:
                union(a, b)
out = pd.DataFrame({"id": list(D)})
out["cluster"] = out.id.map(find)
out["cluster_size"] = out.groupby("cluster").id.transform("size")
out = out[out.id.isin(set(sel.id))]
out.to_csv(os.path.join(R, "derived", "asn_person_clusters.tsv"), sep="\t", index=False)
print("selected applications", len(out), "clusters", out.cluster.nunique(),
      "applications in multi-application clusters", int((out.cluster_size > 1).sum()))

# ======================================================================
# asn_link_deposits.py
#
# Purpose  : Link the two public deposits of the ASN 2016-18 round. The indicator and outcome table of
#            sosgang/coverage_asn is joined to the anonymized DOI lists of Bologna et al. (Zenodo 4603624)
#            on the key (term, level, sector, articles, citations, h-index), keeping only keys that are
#            unique in both tables. Bibliometric sectors only.
# Reads    : data/asn/github_sosgang/indicatoriASN1618.tsv
#              (https://github.com/sosgang/coverage_asn, execution/execution_scripts/indicatoriASN1618.tsv)
#            data/asn/zenodo_4603624/input-anonymized/ND_INPUT_DATA_ANONYMIZED.tsv
#              (https://doi.org/10.5281/zenodo.4603624, input-anonymized.zip)
# Produces : data/asn/derived/asn1618_bibl_linked_noNames.tsv (one row per linked application; no names)
# Run      : from the repository root -- python R/asn_link_deposits.py
# Notes    : The indicator table contains candidate names. The output keeps the CV identifier, which links
#            back to them, so neither file should be redistributed.
# ======================================================================
import csv, os
import pandas as pd

csv.field_size_limit(10**9)
ROOT = os.path.join("data", "asn")
OUT = os.path.join(ROOT, "derived")
os.makedirs(OUT, exist_ok=True)

ind = pd.read_csv(os.path.join(ROOT, "github_sosgang", "indicatoriASN1618.tsv"), sep="\t", dtype=str)
ib = ind[ind.bibl == "B"].copy().fillna("")

path = os.path.join(ROOT, "zenodo_4603624", "input-anonymized", "ND_INPUT_DATA_ANONYMIZED.tsv")
with open(path, encoding="utf-8") as fh:
    rows = list(csv.reader(fh, delimiter="\t"))[1:]
bo = pd.DataFrame({
    "quadrimestre": [r[0] for r in rows], "fascia": [r[1] for r in rows], "settore": [r[2] for r in rows],
    "n_doi": [r[4] for r in rows], "dois": [r[5] for r in rows], "dois_invalid": [r[6] for r in rows],
    "I1": [r[7] for r in rows], "I2": [r[8] for r in rows], "I3": [r[9] for r in rows],
    "SC": [r[10] for r in rows], "SSD_flag": [r[11] for r in rows],
    "S1_art": [r[12] for r in rows], "S2_cit": [r[13] for r in rows], "S3_h": [r[14] for r in rows]})

key = ["quadrimestre", "fascia", "settore", "I1", "I2", "I3"]
bo["k"] = bo[key].agg("|".join, axis=1)
ib["k"] = ib[key].agg("|".join, axis=1)
bo_u = bo[bo.k.map(bo.k.value_counts()) == 1]
ib_u = ib[ib.k.map(ib.k.value_counts()) == 1]
m = ib_u.merge(bo_u.drop(columns=key), on="k", how="inner")
print("bibliometric applications", len(ib), "| linked on a unique key", len(m))

cols = ["id", "quadrimestre", "fascia", "settore", "ssd", "I1", "I2", "I3", "idSoglia", "abilitato", "SC",
        "SSD_flag", "S1_art", "S2_cit", "S3_h", "n_doi", "dois"]
m[cols].to_csv(os.path.join(OUT, "asn1618_bibl_linked_noNames.tsv"), sep="\t", index=False)

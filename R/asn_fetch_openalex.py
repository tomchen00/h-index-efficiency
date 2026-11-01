# ======================================================================
# asn_fetch_openalex.py
#
# Purpose  : Retrieve OpenAlex work records (citation counts by year) for the
#            CV DOIs of the ASN 2016-18 applications that are h-decisive
#            (exactly one of the article and citation thresholds is met) and
#            near the h threshold (|I3 - S3| <= 2).
# Reads    : data/asn/derived/asn1618_bibl_linked_noNames.tsv
# Produces : data/asn/oa_raw/pages.jsonl   one line per API response (query, time, headers, results)
#            data/asn/oa_raw/single.jsonl  one line per single-record lookup (DOIs unusable in OR filters, and misses)
#            data/asn/oa_raw/selection.tsv the selected applications
#            data/asn/oa_raw/fetch_log.txt
# Run      : from the repository root -- python R/asn_fetch_openalex.py
#            Reads OPENALEX_API_KEY from the environment if set (key-less quota is 1000 credits/day).
# Notes    : Filtered list calls cost 1 credit; single-record lookups are free. The script resumes from
#            pages.jsonl / single.jsonl, so an interrupted run can be restarted.
# ======================================================================
import ast, json, os, time, datetime, urllib.parse, urllib.request, urllib.error
import pandas as pd

ROOT = os.path.join("data", "asn")
OUT = os.path.join(ROOT, "oa_raw")
os.makedirs(OUT, exist_ok=True)
KEY = os.environ.get("OPENALEX_API_KEY", "")
FIELDS = "id,doi,publication_year,publication_date,type,cited_by_count,counts_by_year"
MAILTO = "anonymous@example.org"  # replace with your own address

L = pd.read_csv(os.path.join(ROOT, "derived", "asn1618_bibl_linked_noNames.tsv"), sep="\t", dtype=str)
num = ["I1", "I2", "I3", "S1_art", "S2_cit", "S3_h", "n_doi"]
L = L.dropna(subset=num)
for c in num:
    L[c] = L[c].astype(float).astype(int)
decisive = (L.I1 >= L.S1_art) ^ (L.I2 >= L.S2_cit)
near = (L.I3 - L.S3_h).abs() <= 2
P = L[decisive & near].copy()
P.drop(columns=["dois"]).to_csv(os.path.join(OUT, "selection.tsv"), sep="\t", index=False)

dois = set()
for s in P.dois:
    for d in ast.literal_eval(s):
        d = d.strip().lower()
        if d.startswith("10."):
            dois.add(d)
batchable = sorted(d for d in dois if "," not in d and "|" not in d)
bset = set(batchable)
single_only = sorted(d for d in dois if d not in bset)

log = open(os.path.join(OUT, "fetch_log.txt"), "a", encoding="utf-8")
def say(msg):
    line = f"{datetime.datetime.utcnow().isoformat()}Z {msg}"
    print(line); log.write(line + "\n"); log.flush()

say(f"applications {len(P)}; unique DOIs {len(dois)}; batchable {len(batchable)}; single-only {len(single_only)}; key {'yes' if KEY else 'no'}")

def get(url):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                hdr = {k: r.headers.get(k) for k in ("X-RateLimit-Limit", "X-RateLimit-Remaining")}
                return json.load(r), hdr
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, {"status": 404}
            say(f"HTTP {e.code} attempt {attempt}"); time.sleep(5 * (attempt + 1))
        except Exception as e:
            say(f"error {e!r} attempt {attempt}"); time.sleep(5 * (attempt + 1))
    return None, {"status": "failed"}

pages_path = os.path.join(OUT, "pages.jsonl")
done = set()
if os.path.exists(pages_path):
    for line in open(pages_path, encoding="utf-8"):
        done.update(json.loads(line)["query_dois"])
todo = [d for d in batchable if d not in done]
say(f"batch DOIs remaining {len(todo)}")
found = set()
with open(pages_path, "a", encoding="utf-8") as fo:
    for i in range(0, len(todo), 50):
        chunk = todo[i:i + 50]
        q = {"filter": "doi:" + "|".join(chunk), "per-page": 200, "select": FIELDS, "mailto": MAILTO}
        if KEY:
            q["api_key"] = KEY
        res, hdr = get("https://api.openalex.org/works?" + urllib.parse.urlencode(q))
        fo.write(json.dumps({"query_dois": chunk, "time": datetime.datetime.utcnow().isoformat() + "Z",
                             "headers": hdr, "results": (res or {}).get("results")}) + "\n")
        if (i // 50) % 100 == 0:
            say(f"batch {i // 50}: remaining credits {hdr.get('X-RateLimit-Remaining')}")
        time.sleep(0.12)

# DOIs not returned by the batch calls, plus DOIs that cannot be OR-filtered: free single-record lookups
for line in open(pages_path, encoding="utf-8"):
    for w in json.loads(line)["results"] or []:
        if w.get("doi"):
            found.add(w["doi"].lower().replace("https://doi.org/", ""))
single_path = os.path.join(OUT, "single.jsonl")
sdone = set()
if os.path.exists(single_path):
    for line in open(single_path, encoding="utf-8"):
        sdone.add(json.loads(line)["doi"])
miss = [d for d in batchable if d not in found] + single_only
miss = [d for d in miss if d not in sdone]
say(f"single lookups {len(miss)}")
with open(single_path, "a", encoding="utf-8") as fo:
    for d in miss:
        q = {"select": FIELDS, "mailto": MAILTO}
        if KEY:
            q["api_key"] = KEY
        url = "https://api.openalex.org/works/" + urllib.parse.quote("https://doi.org/" + d, safe=":/") + "?" + urllib.parse.urlencode(q)
        res, hdr = get(url)
        fo.write(json.dumps({"doi": d, "time": datetime.datetime.utcnow().isoformat() + "Z", "result": res}) + "\n")
        time.sleep(0.12)
say("DONE")

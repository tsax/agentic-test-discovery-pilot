"""Classify deposit availability for the 40 kept papers and download fetchable data files.

Usage (from repo root):  python3 -I deposits/fetch.py [--no-download]

Outputs
  deposits/availability.csv   one row per paper
  deposits/data/<paper_id>/   downloaded files (gitignored); zips are unpacked (data files only)
  deposits/data_manifest.csv  one row per downloaded/extracted data file

Sources
  (a) PLOS supplementary files: every <supplementary-material> in corpus/raw/<id>.xml.
      href pattern in the XML is "info:doi/10.1371/journal.pone.NNNNNNN.s00N"; PLOS serves it at
      https://journals.plos.org/plosone/article/file?id=<doi>.s00N&type=supplementary
      (302 to a presigned storage URL; HEAD fails there, GET works).
  (b) Repository links in the data-availability statement. Those were resolved by hand on
      2026-11-10 through each repository's public API (OSF v2, Dataverse native API, Zenodo, Dryad v2,
      GitHub tree, Kaggle public download); results are recorded in REPO below with the evidence.

Stdlib only. Polite: <= 4 concurrent downloads, files > 50 MB skipped.
"""
import concurrent.futures as cf
import csv
import io
import json
import os
import sys
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEP = os.path.join(ROOT, "deposits")
DATA = os.path.join(DEP, "data")
XL = "{http://www.w3.org/1999/xlink}href"
UA = "amazon-pilot-research/0.1 (deposit availability check)"  # Dataverse 403s UAs containing "python-urllib"
MAX_BYTES = 50 * 1024 * 1024
DATA_EXT = (".csv", ".tsv", ".tab", ".txt", ".xlsx", ".xls", ".dta", ".sav", ".rds", ".rdata", ".rda")

# Repository deposits named in data-availability statements, resolved by hand (see module doc).
# status: fetchable | not_fetchable ; files: (filename, url)
REPO = {
    "10.1371_journal.pone.0155614": dict(
        repo="Dryad doi:10.5061/dryad.9rj26", status="not_fetchable",
        detail="Dryad lists PHE&RoleStress_Final_occupation.Final.sav (530 kB) but the v2 API file "
               "download returned HTTP 401 to an unauthenticated client; not retried with other clients",
        files=[]),
    "10.1371_journal.pone.0219606": dict(
        repo="GitHub hx2A/predictit_plos_one", status="not_fetchable",
        detail="repo holds code only; data/DATA.md says the PredictIt trade/order files must be "
               "obtained from PredictIt", files=[]),
    "10.1371_journal.pone.0250664": dict(
        repo="Kaggle elinemoens/poned2039500-dataset", status="fetchable",
        detail="public Kaggle download endpoint returned a zip without login",
        files=[("kaggle_dataset.zip",
                "https://www.kaggle.com/api/v1/datasets/download/elinemoens/poned2039500-dataset")]),
    "10.1371_journal.pone.0253430": dict(
        repo="OSF wfsv2 (view-only link)", status="fetchable", detail="OSF API lists data.csv",
        files=[("data.csv", "https://api.osf.io/v2/nodes/wfsv2/files/osfstorage/?view_only=db700c5432b743ee9d69ffeefb9209a3#data.csv")]),
    "10.1371_journal.pone.0259963": dict(
        repo="OSF btruw (view-only link)", status="fetchable", detail="OSF API lists 6 CSVs under Data/",
        files="osf:btruw:c4c137fce44548e9a6021e9eb8e9dff4:Data"),
    "10.1371_journal.pone.0273977": dict(
        repo="Harvard Dataverse doi:10.7910/DVN/9CSTGE", status="fetchable",
        detail="fivecanadas_main_data (Stata 14 original, unrestricted)",
        files=[("fivecanadas_main_data.dta",
                "https://dataverse.harvard.edu/api/access/datafile/6429075?format=original")]),
    "10.1371_journal.pone.0281785": dict(
        repo="OSF mjze3", status="fetchable", detail="OSF API lists Data/epistemic_data.xlsx",
        files="osf:mjze3::Data"),
    "10.1371_journal.pone.0289312": dict(
        repo="UK Data Service SN 6614 (+ OSF 3yn2b code)", status="not_fetchable",
        detail="Understanding Society requires UKDS registration/licence; OSF 3yn2b has only Stata .do files",
        files=[]),
    "10.1371_journal.pone.0327786": dict(
        repo="Zenodo doi:10.5281/zenodo.18180981", status="not_fetchable",
        detail="DOI does not resolve (doi.org 404; Zenodo API 404); PLOS S1 Data is used instead",
        files=[]),
    "10.1371_journal.pone.0342063": dict(
        repo="OSF s2knw (view-only link)", status="fetchable", detail="OSF API lists VacClean.xlsx",
        files=[("VacClean.xlsx", "https://api.osf.io/v2/nodes/s2knw/files/osfstorage/?view_only=4de7d7296d0a44acab11174e357776b#VacClean.xlsx")]),
    "10.1371_journal.pone.0342588": dict(
        repo="openICPSR project 245012", status="not_fetchable",
        detail="openICPSR download requires login (and returned a Cloudflare challenge, HTTP 403)",
        files=[]),
}

# Papers with no fetchable deposit: classification of the data-availability statement (by hand).
NO_DEPOSIT = {
    "10.1371_journal.pone.0178359": ("not_fetchable", "on request to author; IRB-restricted"),
    "10.1371_journal.pone.0190607": ("not_fetchable", "ethical restriction; de-identified data on request"),
    "10.1371_journal.pone.0215957": ("not_fetchable", "on request via data access committee"),
    "10.1371_journal.pone.0280708": ("not_fetchable", "VA IRB restriction; on request"),
    "10.1371_journal.pone.0316899": ("not_fetchable", "CRRS survey; institutional application required"),
    "10.1371_journal.pone.0320807": ("not_fetchable", "confidential; access via Hama Study office"),
    "10.1371_journal.pone.0335670": ("not_fetchable", "confidential administrative data; Ministry of Finance approval"),
    "10.1371_journal.pone.0272792": ("not_fetchable", "third-party subscription (COMPUSTAT/WRDS)"),
    "10.1371_journal.pone.0283041": ("not_fetchable", "third-party, paid registration (BPS Indonesia)"),
    "10.1371_journal.pone.0284692": ("not_fetchable", "third-party subscription (CSMAR, Wind)"),
    "10.1371_journal.pone.0291310": ("none", "third-party public source (World Development Indicators); no analysis file deposited"),
    "10.1371_journal.pone.0347862": ("none", "third-party source (OECD TALIS 2018); no analysis file deposited"),
    "10.1371_journal.pone.0272519": ("none", "'All relevant data are within the paper'; no SI data file"),
    "10.1371_journal.pone.0282251": ("none", "'All relevant data are within the paper'; no SI data file"),
    "10.1371_journal.pone.0283571": ("none", "statement says data in SI, but SI is a questionnaire (DOCX) and test output (DOC)"),
    "10.1371_journal.pone.0298455": ("none", "statement says data in SI, but SI is a DOCX and a Stata .do file; sources are paid yearbooks/CSMAR"),
}


def plos_supp(pid):
    """Return list of (sid, mimetype, label, caption, url) for supplementary files."""
    root = ET.parse(os.path.join(ROOT, "corpus", "raw", pid + ".xml")).getroot()
    out = []
    for sm in root.iter("supplementary-material"):
        href = sm.get(XL) or ""
        doi = href.replace("info:doi/", "")
        lab = (sm.findtext("label") or "").strip()
        cap = " ".join("".join(sm.find("caption").itertext()).split()) if sm.find("caption") is not None else ""
        url = f"https://journals.plos.org/plosone/article/file?id={doi}&type=supplementary"
        out.append((doi.rsplit(".", 1)[-1], sm.get("mimetype", ""), lab, cap, url))
    return out


def supp_is_data_candidate(mime, lab, cap):
    if lab.lower().startswith("attachment"):  # reviewer-response attachments
        return False
    m = mime.lower()
    return ("spreadsheet" in m or "ms-excel" in m or m in ("application/zip", "text/csv")
            or "(CSV)" in cap or "(DTA)" in cap or "(SAV)" in cap)


EXT_BY_MIME = {"spreadsheetml": ".xlsx", "ms-excel": ".xls", "zip": ".zip", "csv": ".csv"}


def get(url, dest):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=120) as r:
        n = r.headers.get("Content-Length")
        if n and int(n) > MAX_BYTES:
            return f"skipped: {int(n)} bytes > 50 MB"
        buf = io.BytesIO()
        while True:
            chunk = r.read(1 << 16)
            if not chunk:
                break
            buf.write(chunk)
            if buf.tell() > MAX_BYTES:
                return "skipped: > 50 MB while streaming"
    with open(dest, "wb") as f:
        f.write(buf.getvalue())
    return "ok"


def osf_list(node, vo, folder):
    """List files (name, download url) in an OSF top-level folder via the public v2 API."""
    q = f"?view_only={vo}" if vo else ""

    def j(u):
        req = urllib.request.Request(u, headers={"User-Agent": UA})
        return json.load(urllib.request.urlopen(req, timeout=60))

    def walk(u):
        res = []
        while u:
            d = j(u if "view_only" in u or not vo else u + ("&" if "?" in u else "?") + f"view_only={vo}")
            for f in d["data"]:
                a = f["attributes"]
                if a["kind"] == "folder":
                    res += walk(f["relationships"]["files"]["links"]["related"]["href"])
                else:
                    res.append((a["name"], f["links"]["download"], a.get("size") or 0))
            u = d["links"].get("next")
        return res

    top = j(f"https://api.osf.io/v2/nodes/{node}/files/osfstorage/{q}")
    for f in top["data"]:
        a = f["attributes"]
        if a["kind"] == "folder" and a["name"] == folder:
            return walk(f["relationships"]["files"]["links"]["related"]["href"])
        if a["kind"] == "file" and not folder:
            pass
    return []


def osf_file(node, vo, name):
    req = urllib.request.Request(f"https://api.osf.io/v2/nodes/{node}/files/osfstorage/?view_only={vo}",
                                 headers={"User-Agent": UA})
    d = json.load(urllib.request.urlopen(req, timeout=60))
    for f in d["data"]:
        if f["attributes"]["name"] == name:
            return f["links"]["download"]
    raise RuntimeError(f"{name} not on OSF {node}")


def safe_unzip(zpath, outdir):
    """Extract data-like members only, flattening paths; nested zips are opened one level deep."""
    got = []
    with zipfile.ZipFile(zpath) as z:
        for info in z.infolist():
            name = info.filename
            base = os.path.basename(name)
            if not base or "__MACOSX" in name or base.startswith("._"):
                continue
            low = base.lower()
            if info.file_size > MAX_BYTES:
                continue
            if low.endswith(".zip"):
                inner = zipfile.ZipFile(io.BytesIO(z.read(info)))
                for ii in inner.infolist():
                    b2 = os.path.basename(ii.filename)
                    if b2 and b2.lower().endswith(DATA_EXT) and "__MACOSX" not in ii.filename and ii.file_size <= MAX_BYTES:
                        dest = os.path.join(outdir, "unz__" + b2)
                        with open(dest, "wb") as f:
                            f.write(inner.read(ii))
                        got.append((dest, name + "!" + ii.filename))
            elif low.endswith(DATA_EXT):
                dest = os.path.join(outdir, "unz__" + base)
                with open(dest, "wb") as f:
                    f.write(z.read(info))
                got.append((dest, name))
    return got


def main(download=True):
    pids = sorted(f[:-5] for f in os.listdir(os.path.join(ROOT, "corpus", "parsed")) if f.endswith(".json"))
    jobs = []  # (pid, source, filename, url)
    rows = {}
    for pid in pids:
        paper = json.load(open(os.path.join(ROOT, "corpus", "parsed", pid + ".json")))
        supp = plos_supp(pid)
        cand = [s for s in supp if supp_is_data_candidate(s[1], s[2], s[3])]
        routes, details = [], []
        for sid, mime, lab, cap, url in cand:
            ext = next((e for k, e in EXT_BY_MIME.items() if k in mime), ".bin")
            jobs.append((pid, f"PLOS {lab or sid}", f"{sid}{ext}", url))
            routes.append(f"PLOS SI {lab} ({sid}{ext})")
        r = REPO.get(pid)
        if r:
            routes.append(r["repo"])
            details.append(r["detail"])
            if r["status"] == "fetchable":
                if isinstance(r["files"], str):
                    _, node, vo, folder = r["files"].split(":")
                    for name, url, size in osf_list(node, vo, folder):
                        if name.lower().endswith(DATA_EXT):
                            jobs.append((pid, r["repo"], name, url))
                else:
                    for name, url in r["files"]:
                        if "api.osf.io" in url:
                            node = url.split("/nodes/")[1].split("/")[0]
                            vo = url.split("view_only=")[1].split("#")[0]
                            url = osf_file(node, vo, url.split("#")[1])
                        jobs.append((pid, r["repo"], name, url))
        rows[pid] = dict(paper_id=pid, doi=paper["doi"], year=paper["year"],
                         n_plos_si=len(supp), plos_si_data_candidates=len(cand),
                         repo_route=r["repo"] if r else "", routes="; ".join(routes),
                         detail="; ".join(details), da_statement=paper["data_availability"][:300])

    results = {}
    if download:
        os.makedirs(DATA, exist_ok=True)

        def run(job):
            pid, src, name, url = job
            d = os.path.join(DATA, pid)
            os.makedirs(d, exist_ok=True)
            dest = os.path.join(d, name)
            try:
                st = get(url, dest)
            except Exception as e:  # network errors are recorded, not fatal
                st = f"error: {type(e).__name__}: {e}"
            return job, st, dest

        with cf.ThreadPoolExecutor(max_workers=4) as ex:
            for job, st, dest in ex.map(run, jobs):
                results.setdefault(job[0], []).append((job, st, dest))

    manifest = []
    for pid, lst in results.items():
        for (p, src, name, url), st, dest in lst:
            if st != "ok":
                manifest.append(dict(paper_id=pid, source=src, file=name, from_archive="", status=st, url=url))
                continue
            if name.lower().endswith(".zip") and zipfile.is_zipfile(dest):
                inner = safe_unzip(dest, os.path.dirname(dest))
                manifest.append(dict(paper_id=pid, source=src, file=name, from_archive="",
                                     status=f"ok (zip, {len(inner)} data members)", url=url))
                for ip, member in inner:
                    manifest.append(dict(paper_id=pid, source=src, file=os.path.basename(ip),
                                         from_archive=name + ":" + member, status="ok", url=url))
            else:
                manifest.append(dict(paper_id=pid, source=src, file=name, from_archive="", status="ok", url=url))

    # classification
    for pid, row in rows.items():
        files = [m for m in manifest if m["paper_id"] == pid and m["status"] == "ok"
                 and m["file"].lower().endswith(DATA_EXT)]
        if files:
            row["category"] = "fetchable"
            row["data_files"] = "; ".join(m["file"] for m in files)
        else:
            row["data_files"] = ""
            if pid in REPO and REPO[pid]["status"] == "not_fetchable":
                row["category"] = "not_fetchable"
            elif pid in NO_DEPOSIT:
                row["category"] = NO_DEPOSIT[pid][0]
                row["detail"] = "; ".join(x for x in (row["detail"], NO_DEPOSIT[pid][1]) if x)
            elif row["plos_si_data_candidates"]:
                row["category"] = "not_fetchable"
                row["detail"] = "; ".join(x for x in (row["detail"], "SI data file listed and downloaded, but it contains no machine-readable data file (e.g. empty zip)") if x)
            else:
                row["category"] = "UNCLASSIFIED"
    cols = ["paper_id", "doi", "year", "category", "routes", "data_files", "detail",
            "n_plos_si", "plos_si_data_candidates", "da_statement"]
    with open(os.path.join(DEP, "availability.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for pid in pids:
            w.writerow(rows[pid])
    if download:
        with open(os.path.join(DEP, "data_manifest.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["paper_id", "source", "file", "from_archive", "status", "url"])
            w.writeheader()
            w.writerows(manifest)
    cnt = {}
    for r in rows.values():
        cnt[r["category"]] = cnt.get(r["category"], 0) + 1
    print("availability:", cnt)


if __name__ == "__main__":
    main(download="--no-download" not in sys.argv)

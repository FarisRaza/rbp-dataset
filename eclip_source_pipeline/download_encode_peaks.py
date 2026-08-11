import requests
import os

OUTDIR = "encode_peaks"
os.makedirs(OUTDIR, exist_ok=True)

# Step 1: get all eCLIP experiments (confirmed working - returns 252)
params = {
    "type": "Experiment",
    "assay_title": "eCLIP",
    "assembly": "GRCh38",
    "internal_tags": "ENCORE",
    "format": "json",
    "limit": "all"
}
r = requests.get("https://www.encodeproject.org/search/", params=params,
                  headers={"accept": "application/json"})
experiments = r.json()["@graph"]
print(f"Found {len(experiments)} eCLIP experiments")

downloaded = 0
skipped_no_match = 0

for exp in experiments:
    target = exp.get("target", {})
    rbp_name = target.get("label", "unknown") if isinstance(target, dict) else "unknown"

    exp_url = "https://www.encodeproject.org" + exp["@id"] + "?format=json"
    exp_data = requests.get(exp_url, headers={"accept": "application/json"}).json()

    files = exp_data.get("files", [])

    # find the combined/preferred_default peaks file
    match = None
    for f in files:
        if (f.get("output_type") == "peaks" and
            f.get("file_format") == "bed" and
            f.get("preferred_default") is True):
            match = f
            break

    if match is None:
        skipped_no_match += 1
        continue

    file_url = "https://www.encodeproject.org" + match["href"]
    fname = os.path.join(OUTDIR, f"{rbp_name}_{match['accession']}.bed.gz")
    if os.path.exists(fname):
        downloaded += 1
        continue

    resp = requests.get(file_url)
    with open(fname, "wb") as out:
        out.write(resp.content)
    downloaded += 1
    print(f"Downloaded {fname}")

print(f"\nDone. {downloaded} files downloaded to {OUTDIR}/")
print(f"{skipped_no_match} experiments had no preferred_default bed peaks file")


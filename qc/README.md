# Exploratory and sanity-check scripts

Each `check_<family>.py` command reads one completed feature sidecar, runs
family-specific invariants, and writes a Markdown report plus PNG figures under
`WORK_DIR/reports/`. The identifier script reads the catalog directly.

Run a single family from the repository root:

```bash
python qc/check_idr.py --work-dir /data/human-proteome-build
python qc/check_eclip.py --work-dir /data/human-proteome-build
python qc/check_opentargets.py --work-dir /data/human-proteome-build
```

Use `--output-dir` to put reports elsewhere. Add `--strict` in automated tests
or release workflows when a failed invariant should produce a nonzero exit
status:

```bash
python qc/check_ptm.py \
  --work-dir /data/human-proteome-build \
  --output-dir /data/release-qc \
  --strict
```

Available commands are:

- `check_identifiers.py`
- `check_cider.py` (whole sequence, per-IDR, and per-domain CIDER)
- `check_idr.py`
- `check_domains.py`
- `check_go.py`
- `check_eclip.py`
- `check_interpro.py`
- `check_ptm.py`
- `check_opentargets.py`
- `check_cdcode.py`
- `check_string.py`
- `check_go_roles.py`
- `check_pslab.py`
- `check_rcsb.py`

The broader command below remains useful when reports for every available
sidecar are wanted in one run:

```bash
python rbp_pipeline/generate_feature_reports.py \
  --work-dir /data/human-proteome-build \
  --features all \
  --skip-unavailable
```


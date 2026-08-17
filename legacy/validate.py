"""Fidelity checks: recompute values for proteins already in the master table.

Every feature family in this package claims to reproduce the logic that built
the existing 17,582 proteins. This script tests that claim the only way that
counts -- by recomputing values for proteins already in the table and comparing
against what is stored there.

Run it after changing any module, and after refreshing any source file. A
regression here means new rows would be subtly incomparable to old ones.

    python validate.py          # base interpreter: CIDER, domains, GO, STRING
    python validate.py --all    # also spawns isoenv/psenv for IDRs and PSLab

Not covered: Open Targets, whose columns are copied through from source files
rather than computed, and the identity columns, which are reconstructed rather
than reproduced (see README).
"""

import ast
import csv
import subprocess
import sys

import paths
import schema

csv.field_size_limit(1 << 30)

# Proteins to test against. Chosen to span the interesting cases: multiple
# repeats of one domain (Q9NQ94, 3x RRM), several IDRs, condensate membership
# (P49588), and no domains at all.
SAMPLE = ["Q9NQ94", "P49588", "Q9NY61", "P61221", "Q8NE71", "P42765"]


def load_sample(accessions):
    """Pull one row per accession out of the master table."""
    wanted, rows = set(accessions), {}
    with open(paths.MASTER_TABLE, newline="", encoding="utf-8",
              errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        index = {c: i for i, c in enumerate(header)}
        for row in reader:
            accession = row[index["uniprot_id"]]
            if accession in wanted and accession not in rows:
                rows[accession] = row
                if len(rows) == len(wanted):
                    break
    return rows, index


def parse(text):
    if text is None or not text.strip():
        return None
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return text


class Report:
    def __init__(self):
        self.passed = self.failed = 0

    def check(self, label, expected, actual, tolerance=None):
        if tolerance is not None and isinstance(expected, float) \
                and isinstance(actual, float):
            same = abs(expected - actual) <= tolerance * max(1.0, abs(expected))
        else:
            same = expected == actual
        if same:
            self.passed += 1
        else:
            self.failed += 1
            print(f"    FAIL {label}")
            print(f"      stored     {str(expected)[:150]}")
            print(f"      recomputed {str(actual)[:150]}")

    def summary(self, title):
        mark = "ok" if not self.failed else "FAILED"
        print(f"  {title}: {self.passed} matched, {self.failed} differed  [{mark}]")
        return self.failed == 0


def check_cider(rows, index):
    import cider

    report = Report()
    for accession, row in rows.items():
        values = cider.whole_sequence(row[index["sequence"]])
        for column in schema.CIDER_COLUMNS:
            stored = parse(row[index[column]])
            actual = values[column]
            if isinstance(stored, float) and isinstance(actual, float):
                report.check(f"{accession}.{column}", stored, actual, tolerance=1e-9)
            else:
                report.check(f"{accession}.{column}", stored, actual)
    return report.summary("CIDER (whole sequence)")


def check_domains_and_go(rows, index):
    from Bio import SwissProt
    import domains as domains_module
    import go as go_module

    remaining = set(rows)
    computed = {}
    with open(paths.UNIPROT_DAT, "r", encoding="utf-8", errors="replace") as fh:
        for record in SwissProt.parse(fh):
            hits = remaining.intersection(record.accessions)
            if not hits:
                continue
            geometry = domains_module.attach_sequences(
                domains_module.domains_for_record(record), record.sequence)
            terms = go_module.go_for_record(record)
            for accession in hits:
                computed[accession] = (geometry, terms)
            remaining -= hits
            if not remaining:
                break

    report = Report()
    for accession, row in rows.items():
        geometry, terms = computed[accession]
        for column in schema.DOMAIN_GEOMETRY_COLUMNS:
            stored = parse(row[index[column]])
            actual = geometry[column]
            if column == "Domains":
                report.check(f"{accession}.{column}",
                             sorted(stored or []), sorted(actual))
            else:
                report.check(f"{accession}.{column}", stored or {}, actual)
        for column in schema.GO_COLUMNS:
            report.check(f"{accession}.{column}",
                         parse(row[index[column]]) or [], terms[column])
    return report.summary("Domains + GO")


def check_string(rows, index):
    import string_ppi

    ensps = {row[index["ENSP_clean"]] for row in rows.values()
             if row[index["ENSP_clean"]]}
    adjacency, _ = string_ppi.scan(paths.STRING_LINKS, ensps)

    report = Report()
    for accession, row in rows.items():
        ensp = row[index["ENSP_clean"]]
        if not ensp:
            continue
        stored = parse(row[index["PPI_ENSP_Partners"]]) or {}
        report.check(f"{accession}.PPI_ENSP_Partners", stored,
                     adjacency.get(ensp, {}))
    return report.summary("STRING partners")


def check_cdcode(rows, index):
    import cdcode

    problems = cdcode.verify_alignment(paths.KAPPEL)
    if problems:
        print(f"    FAIL CD-CODE file ordering: {len(problems)} member tables "
              f"disagree with their declared protein count")
        print(f"      first few: {problems[:5]}")
        return False

    protein_to_condensates, attributes = cdcode.build_index(paths.KAPPEL)
    report = Report()
    skipped = 0
    for accession, row in rows.items():
        actual = cdcode.lookup(accession, protein_to_condensates, attributes)
        for column in schema.CDCODE_COLUMNS:
            stored = parse(row[index[column]])
            if isinstance(stored, list):
                report.check(f"{accession}.{column}", stored, actual[column])
            else:
                # The numeric CD-CODE columns were written as numpy reprs
                # (``[np.int64(9606)]``) which literal_eval refuses, so there is
                # nothing to compare against. Counted as skipped, not passed.
                skipped += 1
    passed = report.summary("CD-CODE")
    if skipped:
        print(f"    ({skipped} cells skipped: stored as numpy reprs, not parseable)")
    return passed


def main(argv):
    rows, index = load_sample(SAMPLE)
    print(f"validating against {len(rows)} proteins already in the master table\n")

    results = [
        check_cider(rows, index),
        check_domains_and_go(rows, index),
        check_string(rows, index),
        check_cdcode(rows, index),
    ]

    if "--all" in argv:
        print("\n  spawning per-environment checks")
        for interpreter, script in ((paths.PYTHON_METAPREDICT, "validate_idr.py"),
                                    (paths.PYTHON_PSLAB, "validate_pslab.py")):
            outcome = subprocess.run([interpreter, "-u", script])
            results.append(outcome.returncode == 0)

    print()
    if all(results):
        print("ALL CHECKS PASSED")
        return 0
    print("SOME CHECKS FAILED -- new rows may not be comparable to existing ones")
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

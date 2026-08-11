"""Master pipeline: rebuild the reviewed human proteome and NCBI isoforms.

Typical use::

    python setup_environment.py
    python rebuild_from_scratch.py all --work-dir D:/human_isoform_rebuild

Stages are resumable.  ``download`` obtains current reviewed human Swiss-Prot,
the NCBI human gene package and (by default) the current Ensembl peptide FASTA.
``catalog`` creates stable canonical and sequence-deduplicated isoform rows.
``features`` invokes one independent ``annotate_<family>.py`` program per
selected family.  ``assemble`` performs keyed left joins into the final table.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import assemble_features
import ensembl_ids
import isoform_catalog
import ncbi_isoforms
import opentargets
import paths
import swissprot_source
from catalog_io import write_rows
from rebuild_schema import BASE_COLUMNS


HERE = os.path.dirname(os.path.abspath(__file__))

# Dependency order is intentional. CIDER consumes IDR/domain sidecars, GO roles
# consumes GO, and PSLab consumes IDRs.
DEFAULT_FEATURES = [
    "idr",
    "domains",
    "go",
    "eclip",
    "interpro",
    "ptm",
    "opentargets",
    "cdcode",
    "string",
    "go_roles",
    "cider",
    "pslab",
]
OPTIONAL_FEATURES = ["rcsb"]
ALL_FEATURES = DEFAULT_FEATURES + OPTIONAL_FEATURES


def layout(work_dir, opentargets_release=opentargets.DEFAULT_RELEASE):
    root = os.path.abspath(work_dir)
    sources = os.path.join(root, "sources")
    opentargets_root = os.path.join(
        sources, "opentargets", opentargets_release
    )
    return {
        "root": root,
        "sources": sources,
        "swissprot": os.path.join(sources, "human_reviewed_swissprot.dat"),
        "ncbi_zip": os.path.join(sources, "ncbi_human_gene.zip"),
        "ncbi_dir": os.path.join(sources, "ncbi_human_gene"),
        "ensembl": os.path.join(sources, "Homo_sapiens.GRCh38.pep.all.fa.gz"),
        "opentargets_root": opentargets_root,
        "opentargets_associations": os.path.join(
            opentargets_root, "association_by_datatype_direct"
        ),
        "opentargets_expression": os.path.join(
            opentargets_root, "expression"
        ),
        "opentargets_diseases": os.path.join(
            opentargets_root, "disease.parquet"
        ),
        "opentargets_release": opentargets_release,
        "catalog": os.path.join(root, "catalog", "human_protein_isoforms.parquet"),
        "catalog_audit": os.path.join(root, "catalog", "catalog.audit.json"),
        "features": os.path.join(root, "features"),
        "final": os.path.join(root, "human_proteome_isoforms_features.parquet"),
        "run_manifest": os.path.join(root, "run_manifest.json"),
    }


def parse_features(value):
    if not value or value == "default":
        return list(DEFAULT_FEATURES)
    if value == "all":
        return list(ALL_FEATURES)
    selected = [item.strip() for item in value.split(",") if item.strip()]
    unknown = sorted(set(selected).difference(ALL_FEATURES))
    if unknown:
        raise ValueError(f"unknown feature families: {unknown}")
    # Restore dependency order regardless of the user's comma order.
    return [family for family in ALL_FEATURES if family in selected]


def download_stage(config, args, selected=None):
    os.makedirs(config["sources"], exist_ok=True)
    swissprot_source.download_human_reviewed(
        config["swissprot"], force=args.force
    )
    datasets_exe = args.datasets_exe
    if not datasets_exe:
        project_copy = os.path.join(paths.KAPPEL, "datasets.exe")
        datasets_exe = project_copy if os.path.exists(project_copy) else None
    ncbi_isoforms.download_human_gene_package(
        config["ncbi_zip"],
        datasets_exe=datasets_exe,
        include_predicted=args.include_predicted,
        include_orphan_refseq=args.include_orphan_refseq,
        force=args.force,
    )
    ncbi_paths = ncbi_isoforms.extract_gene_package(
        config["ncbi_zip"], config["ncbi_dir"], force=args.force
    )
    if not args.no_ensembl_fallback:
        ensembl_ids.download_current_peptides(config["ensembl"], force=args.force)
    if selected is None or "opentargets" in selected:
        print(f"downloading Open Targets release {args.opentargets_release}")
        opentargets.download_release(
            config["opentargets_root"],
            release=args.opentargets_release,
            force=args.force,
        )
    return ncbi_paths


def _ncbi_paths(config):
    return ncbi_isoforms.extract_gene_package(
        config["ncbi_zip"], config["ncbi_dir"], force=False
    )


def catalog_stage(config, args):
    if os.path.exists(config["catalog"]) and not args.force:
        print(f"catalog exists, keeping: {config['catalog']}")
        return config["catalog"]
    ncbi_paths = _ncbi_paths(config)
    ensembl_index = None
    ensembl_release = None
    if not args.no_ensembl_fallback and os.path.exists(config["ensembl"]):
        print("indexing current Ensembl peptide sequences")
        ensembl_index = ensembl_ids.build_sequence_index(config["ensembl"])
        manifest_path = config["ensembl"] + ".manifest.json"
        if os.path.exists(manifest_path):
            with open(manifest_path, encoding="utf-8") as handle:
                manifest = json.load(handle)
            ensembl_release = manifest.get("resolved_url") or manifest.get("url")
    print("building canonical + sequence-unique NCBI isoform catalog")
    rows, audit = isoform_catalog.build_catalog(
        config["swissprot"],
        ncbi_paths["gene_report"],
        ncbi_paths["product_report"],
        ncbi_paths["protein_fasta"],
        include_predicted=args.include_predicted,
        ensembl_index=ensembl_index,
        ensembl_release=ensembl_release,
    )
    write_rows(rows, config["catalog"], BASE_COLUMNS)
    os.makedirs(os.path.dirname(config["catalog_audit"]), exist_ok=True)
    with open(config["catalog_audit"], "w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2)
    print(json.dumps(audit, indent=2))
    return config["catalog"]


def _feature_sources(config):
    def prefer(downloaded, historical):
        return downloaded if os.path.exists(downloaded) else historical

    return {
        "swissprot": config["swissprot"],
        "eclip_table": os.path.join(paths.KAPPEL, "rbp_master_with_eclip.csv"),
        "interpro_tsv": os.path.join(paths.KAPPEL, "df_np_unique.interpro.tsv"),
        "ptm_csv": os.path.join(paths.KAPPEL, "df_ptm.csv"),
        "opentargets_associations": prefer(
            config["opentargets_associations"], paths.OT_ASSOCIATIONS
        ),
        "opentargets_expression": prefer(
            config["opentargets_expression"], paths.OT_EXPRESSION
        ),
        "opentargets_diseases": prefer(
            config["opentargets_diseases"], paths.OT_DISEASES
        ),
        "cdcode_root": paths.KAPPEL,
        "string_links": paths.STRING_LINKS,
        "pspred_repo": paths.PSPRED_REPO,
        "rcsb_summary": paths.RCSB_SUMMARY,
    }


SOURCE_REQUIREMENTS = {
    "domains": ["swissprot"],
    "go": ["swissprot"],
    "eclip": ["eclip_table"],
    "interpro": ["interpro_tsv"],
    "ptm": ["ptm_csv"],
    "opentargets": [
        "opentargets_associations", "opentargets_expression", "opentargets_diseases"
    ],
    "cdcode": ["cdcode_root"],
    "string": ["string_links"],
    "pslab": ["pspred_repo"],
    "rcsb": ["rcsb_summary"],
}


def _python_for(family):
    candidate = {
        "idr": paths.PYTHON_METAPREDICT,
        "pslab": paths.PYTHON_PSLAB,
    }.get(family)
    return candidate if candidate and os.path.exists(candidate) else sys.executable


def _feature_command(family, config, sources):
    output = os.path.join(config["features"], family + ".parquet")
    command = [
        _python_for(family),
        os.path.join(HERE, f"annotate_{family}.py"),
        "--input", config["catalog"],
        "--output", output,
    ]
    flags = {
        "swissprot": "--swissprot",
        "eclip_table": "--eclip-table",
        "interpro_tsv": "--interpro-tsv",
        "ptm_csv": "--ptm-csv",
        "opentargets_associations": "--opentargets-associations",
        "opentargets_expression": "--opentargets-expression",
        "opentargets_diseases": "--opentargets-diseases",
        "cdcode_root": "--cdcode-root",
        "string_links": "--string-links",
        "pspred_repo": "--pspred-repo",
        "rcsb_summary": "--rcsb-summary",
    }
    for name, flag in flags.items():
        if sources.get(name):
            command.extend([flag, sources[name]])
    dependency_paths = {
        "idr_features": os.path.join(config["features"], "idr.parquet"),
        "domain_features": os.path.join(config["features"], "domains.parquet"),
        "go_features": os.path.join(config["features"], "go.parquet"),
    }
    if family in {"cider", "pslab"}:
        command.extend(["--idr-features", dependency_paths["idr_features"]])
    if family == "cider":
        command.extend(["--domain-features", dependency_paths["domain_features"]])
    if family == "go_roles":
        command.extend(["--go-features", dependency_paths["go_features"]])
    if family == "opentargets":
        command.extend([
            "--opentargets-release", config["opentargets_release"],
            "--opentargets-cache",
            os.path.join(config["features"], "opentargets_source_index.sqlite"),
            "--opentargets-expression-long-output",
            os.path.join(config["features"], "opentargets_expression_long.parquet"),
            "--opentargets-disease-long-output",
            os.path.join(config["features"], "opentargets_disease_long.parquet"),
        ])
    return command, output


def features_stage(config, args, selected):
    os.makedirs(config["features"], exist_ok=True)
    sources = _feature_sources(config)
    completed = []
    for family in selected:
        missing = [
            name for name in SOURCE_REQUIREMENTS.get(family, [])
            if not os.path.exists(sources[name])
        ]
        if missing:
            message = f"{family}: missing source(s) {missing}"
            if args.skip_unavailable:
                print("SKIP " + message)
                continue
            raise FileNotFoundError(message)
        command, output = _feature_command(family, config, sources)
        if os.path.exists(output) and not args.force:
            print(f"{family}: sidecar exists, keeping {output}")
            completed.append(output)
            continue
        print(f"\n[{family}] {' '.join(command[:2])}")
        subprocess.run(command, check=True, cwd=HERE)
        completed.append(output)
    return completed


def assemble_stage(config, args, selected):
    sidecars = [
        os.path.join(config["features"], family + ".parquet")
        for family in selected
        if os.path.exists(os.path.join(config["features"], family + ".parquet"))
    ]
    if not sidecars:
        raise FileNotFoundError("no feature sidecars are available to assemble")
    manifest = assemble_features.assemble(
        config["catalog"], sidecars, args.output or config["final"]
    )
    print(json.dumps(manifest, indent=2))
    return manifest


def write_run_manifest(config, args, selected):
    manifest = {
        "command": args.command,
        "work_dir": config["root"],
        "catalog": config["catalog"],
        "final": args.output or config["final"],
        "features_requested": selected,
        "include_predicted_refseq": args.include_predicted,
        "include_orphan_refseq": args.include_orphan_refseq,
        "ensembl_exact_sequence_fallback": not args.no_ensembl_fallback,
        "opentargets_release": args.opentargets_release,
        "feature_sidecars": {
            family: os.path.join(config["features"], family + ".parquet")
            for family in selected
        },
    }
    os.makedirs(config["root"], exist_ok=True)
    with open(config["run_manifest"], "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)


def parser():
    default_work = os.environ.get(
        "KAPPEL_REBUILD_DIR", os.path.join(paths.KAPPEL, "rebuild_output")
    )
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "command", choices=["download", "catalog", "features", "assemble", "all"]
    )
    p.add_argument("--work-dir", default=default_work)
    p.add_argument(
        "--features", default="default",
        help="default, all, or comma-separated families; add rcsb with --features all",
    )
    p.add_argument("--output", help="final .parquet or .csv path")
    p.add_argument("--datasets-exe", help="path to NCBI datasets executable")
    p.add_argument(
        "--opentargets-release",
        default=opentargets.DEFAULT_RELEASE,
        help="pinned Open Targets Platform release (default: %(default)s)",
    )
    p.add_argument("--include-predicted", action="store_true", help="include XP_ products")
    p.add_argument(
        "--include-orphan-refseq", action="store_true",
        help="include NCBI proteins whose gene has no reviewed Swiss-Prot parent",
    )
    p.add_argument(
        "--no-ensembl-fallback", action="store_true",
        help="do not exact-match sequences against current Ensembl peptide FASTA",
    )
    p.add_argument("--skip-unavailable", action="store_true")
    p.add_argument("--force", action="store_true")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    selected = parse_features(args.features)
    config = layout(args.work_dir, args.opentargets_release)
    write_run_manifest(config, args, selected)

    if args.command in {"download", "all"}:
        download_stage(config, args, selected)
    if args.command in {"catalog", "all"}:
        catalog_stage(config, args)
    if args.command in {"features", "all"}:
        features_stage(config, args, selected)
    if args.command in {"assemble", "all"}:
        assemble_stage(config, args, selected)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

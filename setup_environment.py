"""Create the two extra interpreters and vendor the PSLab predictor.

Most of the pipeline runs on whatever Python you already have (it needs only
localcider, pandas, numpy, biopython, openpyxl). Two families do not:

  * IDRs need **metapredict**, which pulls in torch;
  * PSLab needs **scikit-learn 1.6 exactly**, because the shipped ``.joblib``
    models were pickled with it and will not unpickle cleanly on a newer one.

Those two requirement sets conflict, so each gets its own virtual environment
under a short path (``%TEMP%/isoenvs`` by default -- see the ENV_ROOT note in
paths.py for why it must not live inside this package). Run this once:

    python setup_environment.py

It is idempotent -- already-satisfied steps are skipped -- so it is safe to
re-run after a partial failure. Use ``--force`` to rebuild from scratch.

The metapredict install is not a plain ``pip install``
-------------------------------------------------------
On Windows two things go wrong, and this script works around both:

1. metapredict publishes no Windows wheel, so pip falls back to the source
   distribution -- whose test fixtures include filenames containing ``|``, a
   character Windows forbids. Extraction fails before the build even starts.
   Fix: unpack the sdist ourselves, skipping those fixture files.

2. The sdist then tries to compile a Cython accelerator
   (``backend.cython.domain_definition``), which needs MSVC. Most machines do
   not have it.
   Fix: drop the extension from the build. metapredict ships a pure-Python
   implementation of the same routine and its own test suite
   (``tests/test_cython_v_python.py``) asserts the two agree, so `idr.py`
   reroutes to it via `install_python_fallback()`. Results are identical; only
   speed differs, and at a few thousand proteins that is seconds.

On Linux and macOS the wheel installs normally and neither workaround is used.
"""

import argparse
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request
import venv

import paths

HERE = os.path.dirname(os.path.abspath(__file__))
# Short path on purpose -- see the ENV_ROOT comment in paths.py. torch nests
# files ~160 characters deep and Windows caps full paths at 260.
ENVS = paths.ENV_ROOT
VENDOR = os.path.join(HERE, "vendor")

METAPREDICT_VERSION = "3.0.2"
SKLEARN_PIN = "scikit-learn==1.6"
PSPRED_URL = "https://github.com/KULL-Centre/_2024_buelow_PSpred.git"

BASE_REQUIREMENTS = [
    "localcider", "pandas", "numpy", "biopython", "openpyxl", "requests",
    "pyarrow", "duckdb",
]
PSLAB_REQUIREMENTS = [SKLEARN_PIN, "MDAnalysis", "numba", "joblib",
                      "biopython", "pandas", "numpy", "pyarrow"]

# Windows forbids these in filenames; the metapredict sdist contains some.
ILLEGAL_IN_FILENAME = re.compile(r'[|<>:"?*]')


def run(command, **kwargs):
    print("    $ " + " ".join(str(c) for c in command[:4]) + " ...", flush=True)
    result = subprocess.run(command, **kwargs)
    if result.returncode != 0:
        raise SystemExit(f"command failed: {' '.join(map(str, command))}")
    return result


def interpreter_for(env_dir):
    binary = "Scripts" if os.name == "nt" else "bin"
    name = "python.exe" if os.name == "nt" else "python"
    return os.path.join(env_dir, binary, name)


def make_venv(env_dir, force):
    if force and os.path.isdir(env_dir):
        shutil.rmtree(env_dir)
    python = interpreter_for(env_dir)
    if os.path.exists(python):
        print(f"  venv already present: {env_dir}")
        return python
    print(f"  creating venv: {env_dir}")
    venv.EnvBuilder(with_pip=True).create(env_dir)
    run([python, "-m", "pip", "install", "-q", "--upgrade", "pip"])
    return python


def installed(python, module):
    try:
        return subprocess.run(
            [python, "-c", f"import {module}"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            timeout=180,
        ).returncode == 0
    except subprocess.TimeoutExpired:
        return False


def fetch_metapredict_source(destination):
    """Download the sdist and unpack it, skipping Windows-illegal filenames."""
    url_json = f"https://pypi.org/pypi/metapredict/{METAPREDICT_VERSION}/json"
    metadata = json.load(urllib.request.urlopen(url_json, timeout=60))
    sdist = next(f["url"] for f in metadata["urls"]
                 if f["filename"].endswith(".tar.gz"))
    print(f"  downloading {sdist}")
    payload = urllib.request.urlopen(sdist, timeout=600).read()

    if os.path.isdir(destination):
        shutil.rmtree(destination)
    os.makedirs(destination, exist_ok=True)

    skipped = 0
    with tarfile.open(fileobj=io.BytesIO(payload)) as archive:
        for member in archive.getmembers():
            if ILLEGAL_IN_FILENAME.search(os.path.basename(member.name)):
                skipped += 1
                continue
            archive.extract(member, path=destination, filter="data")
    if skipped:
        print(f"  skipped {skipped} test fixture(s) with Windows-illegal names")

    return os.path.join(destination, f"metapredict-{METAPREDICT_VERSION}")


def disable_cython_extension(source_dir):
    """Replace setup.py with one that does not build the C extension."""
    setup = os.path.join(source_dir, "setup.py")
    with open(setup, "w", encoding="utf-8") as fh:
        fh.write(
            "from setuptools import setup, find_packages\n\n"
            "# Cython extension deliberately omitted: it needs MSVC, which is not\n"
            "# assumed here. metapredict falls back to\n"
            "# backend/domain_definition.py:__build_domains_from_values, which its\n"
            "# own tests assert is identical to the compiled version. See\n"
            "# idr.install_python_fallback().\n"
            "setup(packages=find_packages(), include_package_data=True)\n"
        )
    print("  disabled the Cython extension in setup.py")


def setup_metapredict(force):
    print("\n[1/3] metapredict environment")
    python = make_venv(os.path.join(ENVS, "mp"), force)

    if installed(python, "metapredict") and not force:
        print("  metapredict already installed")
    else:
        # CPU-only torch: this pipeline never uses a GPU and the CUDA wheels are
        # several GB.
        run([python, "-m", "pip", "install", "-q",
             "--index-url", "https://download.pytorch.org/whl/cpu", "torch"])

        try:
            run([python, "-m", "pip", "install", "-q",
                 f"metapredict=={METAPREDICT_VERSION}"])
            print("  metapredict installed from the published distribution")
        except SystemExit:
            print("  standard install failed; falling back to a patched sdist")
            source = fetch_metapredict_source(os.path.join(VENDOR, "metapredict_src"))
            disable_cython_extension(source)
            run([python, "-m", "pip", "install", "-q", source])

    run([python, "-m", "pip", "install", "-q", "localcider", "numpy", "pandas",
         "pyarrow"])
    return python


def setup_pslab(force):
    print("\n[2/3] PSLab environment")
    python = make_venv(os.path.join(ENVS, "ps"), force)
    if installed(python, "sklearn") and installed(python, "MDAnalysis") and not force:
        print("  requirements already installed")
    else:
        run([python, "-m", "pip", "install", "-q"] + PSLAB_REQUIREMENTS)
    # Feature sidecars are Parquet even when the scientific model environment
    # already existed before the clean rebuild pipeline added this dependency.
    if not installed(python, "pyarrow"):
        run([python, "-m", "pip", "install", "-q", "pyarrow"])
    return python


def setup_pspred(force):
    print("\n[3/3] PSLab predictor (models + code)")
    repo = paths.PSPRED_REPO
    if force and os.path.isdir(repo):
        shutil.rmtree(repo)
    if os.path.isdir(os.path.join(repo, "models")):
        print(f"  already present: {repo}")
    else:
        os.makedirs(os.path.dirname(repo), exist_ok=True)
        run(["git", "clone", "--depth", "1", PSPRED_URL, repo])

    # sequence.py and predictor.py import each other by bare module name and
    # read these two by relative path, so they must sit alongside them.
    scripts = os.path.join(repo, "scripts_colab")
    for source, name in (
        (os.path.join(repo, "data", "residues.csv"), "residues.csv"),
        (os.path.join(repo, "models", "svr_model_nu.joblib"), "svr_model_nu.joblib"),
    ):
        target = os.path.join(scripts, name)
        if not os.path.exists(target):
            shutil.copy(source, target)
            print(f"  staged {name} into scripts_colab/")
    return repo


def verify(metapredict_python, pslab_python):
    print("\nverifying")
    ok = True

    for module in BASE_REQUIREMENTS:
        importable = module.replace("biopython", "Bio").replace("-", "_")
        if not installed(sys.executable, importable):
            print(f"  [ ] base interpreter is missing {module}"
                  f"   -- pip install {module}")
            ok = False
    if ok:
        print("  [x] base interpreter has localcider, pandas, numpy, biopython, openpyxl")

    try:
        check = subprocess.run(
            [metapredict_python, "-c",
             "import metapredict, localcider; print(metapredict.__version__)"],
            capture_output=True, text=True, timeout=180,
        )
    except subprocess.TimeoutExpired:
        check = subprocess.CompletedProcess([], 124, "", "import timed out")
    if check.returncode == 0:
        print(f"  [x] metapredict {check.stdout.strip().splitlines()[-1]}")
    else:
        print("  [ ] metapredict environment is not usable")
        ok = False

    try:
        check = subprocess.run(
            [pslab_python, "-c", "import sklearn, MDAnalysis; print(sklearn.__version__)"],
            capture_output=True, text=True, timeout=180,
        )
    except subprocess.TimeoutExpired:
        check = subprocess.CompletedProcess([], 124, "", "import timed out")
    if check.returncode == 0:
        print(f"  [x] scikit-learn {check.stdout.strip().splitlines()[-1]}")
    else:
        print("  [ ] PSLab environment is not usable")
        ok = False

    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true",
                        help="rebuild environments from scratch")
    args = parser.parse_args()

    os.makedirs(ENVS, exist_ok=True)
    os.makedirs(VENDOR, exist_ok=True)

    metapredict_python = setup_metapredict(args.force)
    pslab_python = setup_pslab(args.force)
    setup_pspred(args.force)

    if not verify(metapredict_python, pslab_python):
        print("\nsetup incomplete")
        return 1

    print("\nsetup complete. Next:")
    print("    python paths.py          # confirm the source files are present")
    print("    python run_all.py --list # see the stages")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

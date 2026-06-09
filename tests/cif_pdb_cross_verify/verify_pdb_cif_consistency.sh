#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA_DIR="$ROOT/tests/cif_pdb_cross_verify"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

python - "$ROOT" "$DATA_DIR" "$TMPDIR" <<'PY'
import os
import re
import sys
import traceback
from pathlib import Path

from propka.run import single


ROOT = Path(sys.argv[1])
DATA_DIR = Path(sys.argv[2])
TMPDIR = Path(sys.argv[3])
PDB_IDS = ["1crn", "1ubq", "2ptc", "4hhb", "1ake", "3sgb"]
TOLERANCE = 0.01


def parse_summary(path):
    rows = []
    at_pka = False
    with open(path, "rt") as handle:
        for line in handle:
            if at_pka:
                if line.startswith("---"):
                    break
                match = re.search(r"\d+\.\d+", line[13:])
                if match is None:
                    continue
                rows.append((line[:13].strip(), float(match.group())))
            elif "model-pKa" in line:
                at_pka = True
    return rows


def run_propka(input_path, workdir):
    workdir.mkdir(parents=True)
    cwd = Path.cwd()
    try:
        os.chdir(workdir)
        single(str(input_path.resolve()))
        return parse_summary(workdir / "{0:s}.pka".format(input_path.stem))
    finally:
        os.chdir(cwd)


def compare_rows(pdb_rows, cif_rows):
    mismatches = []
    if len(pdb_rows) != len(cif_rows):
        mismatches.append(("length", len(pdb_rows), len(cif_rows), None))
    for index, (pdb_row, cif_row) in enumerate(zip(pdb_rows, cif_rows)):
        same_label = pdb_row[0] == cif_row[0]
        same_value = abs(pdb_row[1] - cif_row[1]) <= TOLERANCE
        if not same_label or not same_value:
            mismatches.append(
                (index, pdb_row, cif_row, round(cif_row[1] - pdb_row[1], 4)))
    return mismatches


overall_ok = True
for pdb_id in PDB_IDS:
    pdb_path = DATA_DIR / "{0:s}.pdb".format(pdb_id)
    cif_path = DATA_DIR / "{0:s}.cif".format(pdb_id)
    print("=== {0:s} ===".format(pdb_id.upper()))
    try:
        pdb_rows = run_propka(pdb_path, TMPDIR / pdb_id / "pdb")
        cif_rows = run_propka(cif_path, TMPDIR / pdb_id / "cif")
    except Exception as err:
        overall_ok = False
        print("FAILED to run {0:s}: {1!s}".format(pdb_id.upper(), err))
        traceback.print_exc(limit=1)
        continue

    mismatches = compare_rows(pdb_rows, cif_rows)
    print("pdb_rows={0:d} cif_rows={1:d}".format(len(pdb_rows), len(cif_rows)))
    if mismatches:
        overall_ok = False
        print("mismatch_count={0:d}".format(len(mismatches)))
        for item in mismatches[:20]:
            print("  {0!r}".format(item))
    else:
        print("match: labels and pKa values identical within {0:.2f}".format(
            TOLERANCE))

if not overall_ok:
    raise SystemExit("PDB/mmCIF cross verification failed")
print("OVERALL PASS")
PY

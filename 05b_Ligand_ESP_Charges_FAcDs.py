"""
05b — QM-derived (Jaguar ESP) partial charges for the MD-ready ligands.

WHY THIS EXISTS
---------------
The Desmond systems take their ligand charges from OPLS4. A fixed-charge force field cannot represent
fluorine's low polarisability well, and the worst case is exactly the chemistry this project is about:
a perfluoroalkyl carboxylate — a hard-charged head on a long, weakly-polarisable-modelled tail. The
error lands on the α-carbon's electrophilicity, which is the single quantity the SN2 depends on.
Nothing else in the pipeline compensates for it: the tier gates measure geometry, and MM-GBSA and the
MD both inherit the same charges.

WHAT IT DOES
------------
For every MD-selected ligand (and the controls), it extracts the ligand from the PREPARED structure —
so the geometry is the one the MD will actually start from — runs a Jaguar single-point at the CFG
level of theory with icfit=1 (charges fitted to the electrostatic potential), and writes:

    <Run>/5_TopN_and_Preparation/4_Ligand_ESP_Charges/
        <ligand>_ESP.in / .out          the Jaguar job, kept for audit
        <ligand>_ESP_charges.csv        atom, element, ESP charge
        <ligand>_ESP.mae                the ligand with r_m_charge1 set to the ESP charges
        00_ESP_Charges_Summary.csv      every ligand, every atom, in one sheet

The .mae is what System Builder consumes: load it, tick 'Use custom charges' → 'Partial charges from
structure'. The charges are NOT silently injected into the force field — a charge set that changes the
physics has to be a visible, deliberate act by the person building the system.

VERIFIED, NOT ASSUMED
---------------------
The keyword is icfit=1, confirmed by running it: Jaguar prints 'Atomic charges from electrostatic
potential' and the charges sum to the formal charge (fluoroacetate: -1.000). A hand-written input with
'MAEFILE:' fails (ERROR 1191) — Jaguar wants the geometry inline in &zmat, which is why the input is
built through schrodinger.application.jaguar.input rather than by string templating.

USAGE
-----
    python3 05b_Ligand_ESP_Charges_FAcDs.py <run_folder> [--controls] [--procs N]

Runs under the system python; the Jaguar work is dispatched through $SCHRODINGER/run.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

# -------------------------------------------------------------------------------
# CFG — the single source of truth (level of theory, basis, paths)
# -------------------------------------------------------------------------------
import importlib.util

_HERE = Path(__file__).resolve().parent
_CFG_PATH = _HERE / "00_01_Project_Config_FAcDs.py"
if not _CFG_PATH.exists():
    print(f"Error: config not found at {_CFG_PATH}")
    sys.exit(1)
_spec = importlib.util.spec_from_file_location("pcfg", str(_CFG_PATH))
_pcfg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_pcfg)
CFG = _pcfg.CFG() if isinstance(_pcfg.CFG, type) else _pcfg.CFG

SCHRODINGER = Path(os.environ.get("SCHRODINGER", "/opt/schrodinger"))
SCHROD_RUN = str(SCHRODINGER / "run")
JAGUAR = str(SCHRODINGER / "jaguar")

# The worker runs inside $SCHRODINGER/run (its own python), so it is written out as a script rather
# than imported: the two interpreters do not share an environment.
_WORKER = r'''
import sys, csv, re
from pathlib import Path
from schrodinger import structure
from schrodinger.application.jaguar.input import JaguarInput

prepared, out_dir, stem, basis, dft, procs = sys.argv[1:7]
out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

_AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE',
       'PRO','SER','THR','TRP','TYR','VAL','HID','HIE','HIP','ASH','HOH','NA','CL','SPC','T3P'}

st = next(structure.StructureReader(prepared))
lig = None
for mol in st.molecule:
    names = {a.pdbres.strip() for a in mol.atom}
    if not (names & _AA) and len(mol.atom) > 2:
        lig = mol.extractStructure()
        break
if lig is None:
    print('NO_LIGAND'); sys.exit(2)

chg = sum(a.formal_charge for a in lig.atom)
lig_mae = out_dir / (stem + '_lig.mae')
lig.write(str(lig_mae))

"""
The charges are fitted on the PREPARED geometry — the structure the MD actually starts from — not on
an idealised gas-phase optimum. A charge set derived from a different conformer is a charge set for a
different molecule.  icfit=1 fits to the electrostatic potential; verified against Jaguar's own output
('Atomic charges from electrostatic potential'), and the fitted charges sum to the formal charge.
"""
ji = JaguarInput(name=stem + '_ESP')
ji.setStructure(lig)
ji.setValues({'basis': basis, 'dftname': dft, 'molchg': int(chg), 'multip': 1, 'icfit': 1})
in_path = out_dir / (stem + '_ESP.in')
ji.saveAs(str(in_path))
print('WROTE_INPUT', in_path, 'charge', chg, 'atoms', len(lig.atom))
'''

_PARSE = r'''
import sys, csv
from pathlib import Path
from schrodinger import structure

out_file, lig_mae, stem, out_dir = sys.argv[1:5]
out_dir = Path(out_dir)
text = Path(out_file).read_text(errors='ignore')

# Jaguar prints the fitted charges as label/charge row pairs under one header. The block is read to
# its own end rather than to a fixed line budget: a larger ligand wraps the table over several pairs
# and a fixed cut silently drops the tail.
lines = text.splitlines()
labels, charges = [], []
capture = False
for i, ln in enumerate(lines):
    if 'Atomic charges from electrostatic potential' in ln:
        capture = True
        continue
    if capture:
        s = ln.strip()
        if s.startswith('Atom') and i + 1 < len(lines) and lines[i + 1].strip().startswith('Charge'):
            labels += ln.split()[1:]
            charges += lines[i + 1].split()[1:]
        elif s and not s.startswith('Charge') and labels:
            break
if not labels:
    print('NO_ESP_BLOCK'); sys.exit(3)

vals = [float(c) for c in charges]
rows = list(zip(labels, vals))

csv_path = out_dir / (stem + '_ESP_charges.csv')
with csv_path.open('w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['atom_label', 'esp_charge'])
    w.writerows(rows)

# The charges are written onto the ligand structure as partial charges, which is the form System
# Builder reads ('Use custom charges' -> 'Partial charges from structure').
st = next(structure.StructureReader(lig_mae))
if len(st.atom) != len(vals):
    print(f'ATOM_COUNT_MISMATCH {len(st.atom)} vs {len(vals)}'); sys.exit(4)
for a, q in zip(st.atom, vals):
    a.partial_charge = float(q)
mae_path = out_dir / (stem + '_ESP.mae')
st.write(str(mae_path))

print(f'OK {stem} atoms={len(vals)} sum={sum(vals):+.4f} -> {mae_path.name}')
'''


def _sh(script: str, *args: str, timeout: int = 900) -> subprocess.CompletedProcess:
    """Run a snippet under $SCHRODINGER/run (its interpreter, not ours)."""
    return subprocess.run([SCHROD_RUN, "python3", "-c", script, *args],
                          capture_output=True, text=True, timeout=timeout)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="QM (Jaguar ESP) partial charges for the MD-ready ligands")
    ap.add_argument("run", help="Boltz-2 run folder")
    ap.add_argument("--controls", action="store_true",
                    help="also charge the six control ligands")
    ap.add_argument("--procs", type=int, default=int(getattr(CFG, "QSITE_PROCS", 4)),
                    help="Jaguar processors per job")
    args = ap.parse_args()

    run = Path(args.run).resolve()
    prep_dir = run / "5_TopN_and_Preparation" / "2_Prepared_PDBs"
    if not prep_dir.exists():
        print(f"Error: prepared structures not found: {prep_dir}")
        return 1

    out_dir = run / "5_TopN_and_Preparation" / "4_Ligand_ESP_Charges"
    out_dir.mkdir(parents=True, exist_ok=True)

    structures = sorted(prep_dir.glob("*_Prepared.pdb"))
    if not args.controls:
        structures = [p for p in structures if "_Control_" not in p.name]
    if not structures:
        print("Error: no prepared structures to charge")
        return 1

    basis = str(getattr(CFG, "ESP_CHARGE_BASIS", "6-31G**"))
    dft = str(getattr(CFG, "ESP_CHARGE_DFT", "b3lyp"))

    print(f"  ESP charges — {len(structures)} ligand(s)  |  {dft}/{basis}  |  {args.procs} procs")
    print(f"  output: {out_dir}")

    summary: list[dict] = []
    ok = 0
    for pdb in structures:
        stem = pdb.name.replace("_Prepared.pdb", "")
        print(f"\n  ── {stem}")

        r = _sh(_WORKER, str(pdb), str(out_dir), stem, basis, dft, str(args.procs))
        if "WROTE_INPUT" not in r.stdout:
            print(f"     ! input build failed: {(r.stdout + r.stderr).strip()[:160]}")
            continue
        print(f"     {r.stdout.strip().splitlines()[-1]}")

        in_path = out_dir / f"{stem}_ESP.in"
        jr = subprocess.run([JAGUAR, "run", "-WAIT", "-NOJOBID", in_path.name],
                            cwd=str(out_dir), capture_output=True, text=True, timeout=7200)
        out_path = out_dir / f"{stem}_ESP.out"
        if jr.returncode != 0 or not out_path.exists():
            print(f"     ! Jaguar failed (rc={jr.returncode}); see {out_path.name}")
            continue

        pr = _sh(_PARSE, str(out_path), str(out_dir / f"{stem}_lig.mae"), stem, str(out_dir))
        line = (pr.stdout or pr.stderr).strip().splitlines()[-1] if (pr.stdout or pr.stderr) else ""
        print(f"     {line}")
        if not line.startswith("OK"):
            continue
        ok += 1

        csv_path = out_dir / f"{stem}_ESP_charges.csv"
        with csv_path.open() as fh:
            for row in csv.DictReader(fh):
                summary.append({"structure": stem,
                                "atom_label": row["atom_label"],
                                "esp_charge": row["esp_charge"]})

    if summary:
        s_path = out_dir / "00_ESP_Charges_Summary.csv"
        with s_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["structure", "atom_label", "esp_charge"])
            w.writeheader()
            w.writerows(summary)
        print(f"\n  ✔ {ok}/{len(structures)} ligands charged → {s_path.name}")

    print("""
  NEXT STEP (System Builder, by hand — deliberately not automated):
    load <stem>_ESP.mae, tick 'Use custom charges' → 'Partial charges from structure'.
  A charge set that changes the physics should be applied by a person who knows they are doing it,
  not injected silently by a script.
""")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

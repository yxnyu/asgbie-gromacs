#!/usr/bin/env bash
# ============================================================
# Step 6: ASGBIE Free Energy Calculation
# Uses gmx_MMPBSA to replace AMBER's commercial MMPBSA.py
#
# ASGBIE = Alanine Scanning (AS) +
#          Generalized Born  (GB) +
#          Interaction Entropy (IE)
#
# This script runs two phases:
#   Phase 1: GB + IE baseline binding free energy
#   Phase 2: Alanine scanning (per-residue loop, with resume support)
#
# Core technology:
#   gmx_MMPBSA internally converts GROMACS topology/trajectory to AMBER format,
#   then calls free AmberTools MMPBSA.py for the calculation
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "===== [Step 6] ASGBIE Free Energy Calculation ====="
cd "${WORKDIR}"
mkdir -p asgbie_calc && cd asgbie_calc

# ---------- 6.1 Prepare index file ----------
echo "[6.1] Creating receptor/ligand index groups ..."

# Auto-detect receptor and ligand groups
python3 << 'PYEOF'
import subprocess, os

rec_chains = os.environ.get("RECEPTOR_CHAINS", "A").split()
lig_chains = os.environ.get("LIGAND_CHAINS", "B").split()

cmds = []
for ch in rec_chains:
    cmds.append(f"chain {ch}")
if len(rec_chains) > 1:
    cmds.append(" | ".join([f"chain {c}" for c in rec_chains]))
    cmds.append(f'name {len(cmds)-1} Receptor')

for ch in lig_chains:
    cmds.append(f"chain {ch}")
if len(lig_chains) > 1:
    cmds.append(" | ".join([f"chain {c}" for c in lig_chains]))
    cmds.append(f'name {len(cmds)-1} Ligand')

cmds.append("q")
cmd_str = "\n".join(cmds)

with open("make_ndx_input.txt", "w") as f:
    f.write(cmd_str)
print(f"  Index commands written to make_ndx_input.txt")
PYEOF

gmx make_ndx -f ../prod.tpr -o index_mmpbsa.ndx < make_ndx_input.txt 2>/dev/null || \
    cp ../index.ndx index_mmpbsa.ndx 2>/dev/null || true

# ======================================================================
# Phase 1: GB + IE Baseline Binding Free Energy (without alanine scanning)
# ======================================================================
echo ""
echo "╔══════════════════════════════════════════╗"
echo "║  Phase 1: GB + IE Binding Free Energy    ║"
echo "╚══════════════════════════════════════════╝"

cat > mmpbsa_gbie.in << MMPBSA_EOF
&general
  sys_name             = "ASGBIE-GROMACS-GBIE",
  startframe           = ${MMPBSA_START_FRAME},
  endframe             = ${MMPBSA_END_FRAME},
  interval             = ${MMPBSA_INTERVAL},
  interaction_entropy  = 1,
  ie_segment           = ${IE_SEGMENT},
  temperature          = ${TEMPERATURE},
/

&gb
  igb    = ${IGB},
  saltcon = 0.15,
/
MMPBSA_EOF

echo "[6.2] Running GB + IE calculation ..."
gmx_MMPBSA -O \
    -i mmpbsa_gbie.in \
    -cs ../prod.tpr \
    -ct ../prod.xtc \
    -ci index_mmpbsa.ndx \
    -cg "${REC_GROUP:-1}" "${LIG_GROUP:-13}" \
    -cp ../topol.top \
    -o FINAL_RESULTS_GBIE.dat \
    -eo FINAL_RESULTS_GBIE.csv \
    -nogui 2>&1 | tee gmx_mmpbsa_gbie.log

echo "[Phase 1] GB + IE calculation done ✓"

# ======================================================================
# Phase 2: Alanine Scanning (per-residue loop)
# ======================================================================
echo ""
echo "╔══════════════════════════════════════════╗"
echo "║  Phase 2: Alanine Scanning               ║"
echo "╚══════════════════════════════════════════╝"

# Skip if alanine scanning is disabled
if [[ "${ALA_SCAN_ENABLED:-true}" == "false" ]]; then
    echo "[SKIP] Alanine scanning disabled (ALA_SCAN_ENABLED=false)"
    echo "[Step 6] Done ✓"
    exit 0
fi

# Call batch alanine scanning script
bash "${SCRIPT_DIR}/06b_ala_scan_batch.sh"

# ---------- Aggregate all results ----------
echo ""
echo "[6.3] Aggregating ASGBIE results ..."

# Merge GB+IE and alanine scanning results
python3 << 'PYEOF'
import os, glob, csv, shutil

results_dir = "ala_scan_results"
summary_file = "FINAL_RESULTS_ASGBIE.csv"

# Collect all per-residue scanning results
ala_results = []
if os.path.isdir(results_dir):
    for f in sorted(glob.glob(os.path.join(results_dir, "*.dat"))):
        resname = os.path.basename(f).replace("ala_", "").replace(".dat", "")
        with open(f) as fh:
            content = fh.read()
        # Parse ddG
        for line in content.split("\n"):
            if "DELTA" in line and "TOTAL" in line:
                parts = line.split()
                try:
                    ddg = float(parts[-1])
                    ala_results.append({"residue": resname, "ddG": ddg, "file": f})
                except (ValueError, IndexError):
                    pass

if ala_results:
    # Sort by ddG descending (hotspot residues first)
    ala_results.sort(key=lambda x: x["ddG"], reverse=True)

    with open("ala_scan_summary.csv", "w", newline="") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=["rank", "residue", "ddG_kcal_mol", "hotspot"])
        writer.writeheader()
        for i, r in enumerate(ala_results, 1):
            writer.writerow({
                "rank": i,
                "residue": r["residue"],
                "ddG_kcal_mol": f"{r['ddG']:.2f}",
                "hotspot": "YES" if r["ddG"] > 1.0 else "no",
            })
    print(f"  Alanine scanning summary: ala_scan_summary.csv ({len(ala_results)} residues)")
    print(f"  Hotspot residues (ddG > 1.0 kcal/mol):")
    for r in ala_results:
        if r["ddG"] > 1.0:
            print(f"    {r['residue']:>10s}  ddG = {r['ddG']:+.2f} kcal/mol  ★")
else:
    print("  [INFO] No alanine scanning result files found")

# Copy GB+IE results as main results
if os.path.exists("FINAL_RESULTS_GBIE.dat"):
    shutil.copy("FINAL_RESULTS_GBIE.dat", "FINAL_RESULTS_ASGBIE.dat")
    shutil.copy("FINAL_RESULTS_GBIE.csv", "FINAL_RESULTS_ASGBIE.csv")

PYEOF

echo ""
echo "[Step 6] Done ✓"
echo "  GB+IE results:            FINAL_RESULTS_GBIE.dat/csv"
echo "  Alanine scanning summary: ala_scan_summary.csv"
echo "  Log: gmx_mmpbsa_gbie.log"

#!/usr/bin/env bash
# ============================================================
# Step 7: Result Collection
# Replaces the original AMBER 6result_asgbie.sh
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "===== [Step 7] Result Collection ====="
cd "${WORKDIR}"
mkdir -p "${RESULTS_DIR}"

# ---------- Copy core results ----------
echo "[7.1] Collecting result files ..."

# GB+IE results
cp asgbie_calc/FINAL_RESULTS_GBIE.dat  "${RESULTS_DIR}/" 2>/dev/null || true
cp asgbie_calc/FINAL_RESULTS_GBIE.csv  "${RESULTS_DIR}/" 2>/dev/null || true
cp asgbie_calc/FINAL_RESULTS_ASGBIE.dat "${RESULTS_DIR}/" 2>/dev/null || true
cp asgbie_calc/FINAL_RESULTS_ASGBIE.csv "${RESULTS_DIR}/" 2>/dev/null || true
cp asgbie_calc/gmx_mmpbsa_gbie.log     "${RESULTS_DIR}/" 2>/dev/null || true

# Alanine scanning summary
cp asgbie_calc/ala_scan_summary.csv     "${RESULTS_DIR}/" 2>/dev/null || true

# Alanine scanning per-residue results (archive)
if [[ -d asgbie_calc/ala_scan_results ]] && ls asgbie_calc/ala_scan_results/*.dat &>/dev/null; then
    mkdir -p "${RESULTS_DIR}/ala_scan_details"
    cp asgbie_calc/ala_scan_results/*.dat "${RESULTS_DIR}/ala_scan_details/" 2>/dev/null || true
    cp asgbie_calc/ala_scan_results/*.csv "${RESULTS_DIR}/ala_scan_details/" 2>/dev/null || true
    echo "  Alanine scanning details: ${RESULTS_DIR}/ala_scan_details/"
fi

# RMSD
cp rmsd_backbone.xvg "${RESULTS_DIR}/" 2>/dev/null || true
cp rmsd_protein.xvg  "${RESULTS_DIR}/" 2>/dev/null || true

# Energy
cp em_potential.xvg  "${RESULTS_DIR}/" 2>/dev/null || true
cp nvt_temp.xvg      "${RESULTS_DIR}/" 2>/dev/null || true
cp npt_density.xvg   "${RESULTS_DIR}/" 2>/dev/null || true

# ---------- Check key files ----------
if [[ ! -f "${RESULTS_DIR}/FINAL_RESULTS_GBIE.dat" ]] && [[ ! -f "${RESULTS_DIR}/FINAL_RESULTS_ASGBIE.dat" ]]; then
    echo "[WARN] GB+IE result files not found — please confirm Step 6 has completed"
fi

# ---------- Parse ASGBIE results ----------
echo "[7.2] Parsing ASGBIE results ..."

python3 << 'PYEOF'
import os

results_dir = os.environ.get("RESULTS_DIR", "asgbie_results")

# Print GB+IE results
for dat_name in ["FINAL_RESULTS_ASGBIE.dat", "FINAL_RESULTS_GBIE.dat"]:
    dat_file = os.path.join(results_dir, dat_name)
    if os.path.exists(dat_file):
        print(f"\n{'='*60}")
        print(f"  ASGBIE Result Summary ({dat_name})")
        print(f"{'='*60}")
        with open(dat_file) as f:
            content = f.read()
        for line in content.split('\n'):
            line_s = line.strip()
            if any(kw in line_s for kw in [
                'DELTA', 'TOTAL', 'Binding', 'IE', 'Entropy',
                'VDWAALS', 'EEL', 'EGB', 'ESURF',
                'Mutant', 'Difference'
            ]):
                print(f"  {line_s}")
        print(f"{'='*60}\n")
        break

# Print alanine scanning results
ala_csv = os.path.join(results_dir, "ala_scan_summary.csv")
if os.path.exists(ala_csv):
    import csv
    print(f"\n{'='*60}")
    print(f"  Alanine Scanning Hotspot Residues")
    print(f"{'='*60}")
    with open(ala_csv) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    hotspots = [r for r in rows if float(r.get("ddG_kcal_mol", 0)) > 1.0]
    print(f"  Total scanned: {len(rows)}")
    print(f"  Hotspots:      {len(hotspots)}")
    if hotspots:
        print(f"\n  {'Rank':>4}  {'Residue':>12}  {'ddG (kcal/mol)':>16}  Hotspot")
        print(f"  {'-'*4}  {'-'*12}  {'-'*16}  {'-'*7}")
        for r in hotspots:
            print(f"  {r['rank']:>4}  {r['residue']:>12}  {float(r['ddG_kcal_mol']):>+16.2f}  ★")
    print(f"{'='*60}\n")

# List CSV result files
for csv_name in ["FINAL_RESULTS_ASGBIE.csv", "FINAL_RESULTS_GBIE.csv"]:
    csv_file = os.path.join(results_dir, csv_name)
    if os.path.exists(csv_file):
        print(f"  Detailed CSV results: {csv_file}")
        break
PYEOF

echo "[Step 7] Done ✓"
echo "  All results collected to: ${RESULTS_DIR}/"
ls -lh "${RESULTS_DIR}/" 2>/dev/null

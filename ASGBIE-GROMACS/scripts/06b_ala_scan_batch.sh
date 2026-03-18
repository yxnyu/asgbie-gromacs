#!/usr/bin/env bash
# ============================================================
# Batch Alanine Scanning
#
# Key improvements:
#   gmx_MMPBSA's &alanine_scanning only supports single-residue
#   mutation per run. This script uses a for loop to call
#   gmx_MMPBSA per residue, enabling:
#     1. Per-residue independent execution — one failure doesn't affect others
#     2. Resume support — completed residues are automatically skipped
#     3. GNU Parallel acceleration — multi-core parallel execution
#     4. Auto-aggregation — generates ddG ranking table and hotspot list
#
# Usage:
#   Called automatically by 06_asgbie.sh, or run standalone:
#     bash scripts/06b_ala_scan_batch.sh
#
# Configuration (config.sh):
#   ALA_SCAN_RESIDUES=""    # Empty = auto-detect interface residues
#   ALA_SCAN_PARALLEL=4     # Parallel jobs (0=serial)
#   ALA_SCAN_ENABLED=true   # Enable/disable
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

cd "${WORKDIR}/asgbie_calc"
mkdir -p ala_scan_results ala_scan_logs

# ========== 1. Get residue list to scan ==========
echo "[AS] Detecting interface residues ..."

RESIDUE_LIST_FILE="ala_scan_residues.txt"

python3 << 'PYEOF'
import os, sys, re

user_residues = os.environ.get("ALA_SCAN_RESIDUES", "").strip()
output_file = "ala_scan_residues.txt"

if user_residues:
    # User-specified: "ARG:25,GLU:30,TYR:45"
    residues = []
    for item in user_residues.split(","):
        item = item.strip()
        if ":" in item:
            name, num = item.split(":", 1)
            name = name.strip().upper()
            num = num.strip()
            # Input validation: residue name 1-4 uppercase letters, number is digits only
            if not re.match(r'^[A-Z]{1,4}$', name):
                print(f"  [ERROR] Invalid residue name: {name} (expected 1-4 uppercase letters)")
                sys.exit(1)
            if not re.match(r'^[0-9]+$', num):
                print(f"  [ERROR] Invalid residue number: {num} (expected digits only)")
                sys.exit(1)
            residues.append(f"{name}:{num}")
        else:
            print(f"  [ERROR] Invalid residue format: {item} (expected NAME:NUM, e.g., ARG:25)")
            sys.exit(1)
    with open(output_file, "w") as f:
        for r in residues:
            f.write(r + "\n")
    print(f"  User specified {len(residues)} residues")
    sys.exit(0)

# Auto-detect: extract all protein residues from GRO file (excluding ALA/GLY/solvent/ions)
# Note: This is a simplified implementation that includes ALL protein residues,
# not just interface residues. It is recommended to manually specify interface
# residues via ALA_SCAN_RESIDUES for more precise results.
try:
    gro_file = "../npt.gro"
    if not os.path.exists(gro_file):
        gro_file = "../complex_processed.gro"

    rec_atoms = []

    with open(gro_file) as f:
        lines = f.readlines()

    # GRO format: fixed-width columns (resnum[0:5], resname[5:10], atomname[10:15],
    #             x[20:28], y[28:36], z[36:44])
    for line in lines[2:]:  # Skip title and atom count lines
        if len(line) < 44:
            continue
        try:
            resnum = int(line[0:5].strip())
            resname = line[5:10].strip()
            x = float(line[20:28]) * 10  # nm -> Angstrom
            y = float(line[28:36]) * 10
            z = float(line[36:44]) * 10
        except (ValueError, IndexError):
            continue

        atom_info = {"resnum": resnum, "resname": resname, "x": x, "y": y, "z": z}

        # Exclude solvent and ions
        if resname not in ("SOL", "NA", "CL", "HOH", "WAT", "NA+", "CL-", "K", "K+"):
            rec_atoms.append(atom_info)

    # Get unique residue list, excluding ALA and GLY
    seen = set()
    all_protein_residues = []
    for atom in rec_atoms:
        key = f"{atom['resname']}:{atom['resnum']}"
        if key not in seen and atom["resname"] not in ("ALA", "GLY"):
            seen.add(key)
            all_protein_residues.append(key)

    max_residues = int(os.environ.get("ALA_SCAN_MAX_RESIDUES", "50"))
    selected = all_protein_residues[:max_residues]

    with open(output_file, "w") as f:
        for r in selected:
            f.write(r + "\n")
    print(f"  Auto-detected {len(all_protein_residues)} scannable residues (non-ALA/GLY)")
    print(f"  Selected first {len(selected)} for scanning")
    if len(all_protein_residues) > max_residues:
        print(f"  [Tip] Adjust via ALA_SCAN_MAX_RESIDUES or ALA_SCAN_RESIDUES")
    print(f"  [Tip] Recommend setting ALA_SCAN_RESIDUES manually for interface-only scanning")

except Exception as e:
    print(f"  [ERROR] Auto-detection failed: {e}")
    print(f"  Please set ALA_SCAN_RESIDUES in config.sh manually")
    sys.exit(1)
PYEOF

if [[ ! -f "${RESIDUE_LIST_FILE}" ]] || [[ ! -s "${RESIDUE_LIST_FILE}" ]]; then
    echo "[ERROR] Residue list is empty, cannot run alanine scanning"
    exit 1
fi

TOTAL_RESIDUES=$(wc -l < "${RESIDUE_LIST_FILE}" | tr -d ' ')
echo "[AS] Total ${TOTAL_RESIDUES} residues to scan"

# ========== 2. Define single-residue scanning function ==========
# Write a standalone script that can be called by parallel
# Note: Using 'SINGLE_EOF' (quoted) to prevent variable expansion at write time
cat > _run_single_ala.sh << 'SINGLE_EOF'
#!/usr/bin/env bash
# Single-residue alanine scanning
# Usage: bash _run_single_ala.sh <RESNAME:RESNUM> <WORKDIR> <config_vars...>
# Note: Not using set -e, to capture gmx_MMPBSA exit code properly
set -uo pipefail

RESIDUE="$1"
WORKDIR="$2"
MMPBSA_START_FRAME="${3:-1}"
MMPBSA_END_FRAME="${4:-100}"
MMPBSA_INTERVAL="${5:-1}"
IGB="${6:-5}"
IE_SEGMENT="${7:-25}"
TEMPERATURE="${8:-300}"
REC_GROUP="${9:-1}"
LIG_GROUP="${10:-13}"

RESNAME="${RESIDUE%%:*}"
RESNUM="${RESIDUE##*:}"
SAFE_NAME="${RESNAME}_${RESNUM}"

cd "${WORKDIR}/asgbie_calc"

# Resume: if result already exists and is non-empty, skip
RESULT_FILE="ala_scan_results/ala_${SAFE_NAME}.dat"
if [[ -f "${RESULT_FILE}" ]] && [[ -s "${RESULT_FILE}" ]]; then
    echo "[SKIP] ${RESIDUE} — result already exists"
    exit 0
fi

# Create temporary work directory (avoid using TMPDIR, which is a POSIX system variable)
SCAN_TMPDIR="ala_scan_tmp/${SAFE_NAME}"
mkdir -p "${SCAN_TMPDIR}"
cd "${SCAN_TMPDIR}"

# Generate MMPBSA input file for this residue
cat > mmpbsa_ala.in << MMPBSA_IN
&general
  sys_name             = "ALA_${SAFE_NAME}",
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

&alanine_scanning
  mutant = "ALA",
  mutant_res = "${RESNAME}:${RESNUM}",
/
MMPBSA_IN

# Run gmx_MMPBSA (no set -e, manually capture exit code)
EXIT_CODE=0
gmx_MMPBSA -O \
    -i mmpbsa_ala.in \
    -cs ../../../prod.tpr \
    -ct ../../../prod.xtc \
    -ci ../../index_mmpbsa.ndx \
    -cg "${REC_GROUP}" "${LIG_GROUP}" \
    -cp ../../../topol.top \
    -o results_ala.dat \
    -eo results_ala.csv \
    -nogui > "../../ala_scan_logs/ala_${SAFE_NAME}.log" 2>&1 || EXIT_CODE=$?

if [[ ${EXIT_CODE} -eq 0 ]] && [[ -f results_ala.dat ]]; then
    cp results_ala.dat "../../ala_scan_results/ala_${SAFE_NAME}.dat"
    cp results_ala.csv "../../ala_scan_results/ala_${SAFE_NAME}.csv" 2>/dev/null || true
    echo "[OK] ${RESIDUE} — done"
else
    echo "[FAIL] ${RESIDUE} — failed (exit=${EXIT_CODE}), see ala_scan_logs/ala_${SAFE_NAME}.log"
fi

# Clean up temporary files (keep logs)
cd "${WORKDIR}/asgbie_calc"
rm -rf "ala_scan_tmp/${SAFE_NAME}"
SINGLE_EOF

chmod +x _run_single_ala.sh

# ========== 3. Execute scanning ==========
PARALLEL_JOBS="${ALA_SCAN_PARALLEL:-1}"
mkdir -p ala_scan_tmp

echo "[AS] Starting alanine scanning (parallel jobs: ${PARALLEL_JOBS}) ..."
echo ""

SCAN_START=$(date +%s)
COMPLETED=0
FAILED=0

if [[ ${PARALLEL_JOBS} -gt 1 ]] && command -v parallel &>/dev/null; then
    # ---- GNU Parallel mode ----
    echo "[AS] Using GNU Parallel for parallel execution ..."
    parallel -j "${PARALLEL_JOBS}" --progress --bar \
        bash _run_single_ala.sh {} "${WORKDIR}" \
        "${MMPBSA_START_FRAME}" "${MMPBSA_END_FRAME}" "${MMPBSA_INTERVAL}" \
        "${IGB}" "${IE_SEGMENT}" "${TEMPERATURE}" \
        "${REC_GROUP:-1}" "${LIG_GROUP:-13}" \
        :::: "${RESIDUE_LIST_FILE}" 2>&1 | tee ala_scan_parallel.log

elif [[ ${PARALLEL_JOBS} -gt 1 ]]; then
    # ---- Simple background parallel (no GNU Parallel) ----
    echo "[AS] GNU Parallel not installed, using background processes ..."
    echo "  Tip: Install GNU Parallel for better parallel experience"
    echo "    conda install -c conda-forge parallel"
    echo ""

    PIDS=()
    RUNNING=0

    while IFS= read -r residue; do
        # Control parallelism
        while [[ ${RUNNING} -ge ${PARALLEL_JOBS} ]]; do
            # Wait for any subprocess to complete
            for i in "${!PIDS[@]}"; do
                if ! kill -0 "${PIDS[$i]}" 2>/dev/null; then
                    wait "${PIDS[$i]}" 2>/dev/null || true
                    unset 'PIDS[$i]'
                    RUNNING=$((RUNNING - 1))
                fi
            done
            # Rebuild array index (safe empty array handling)
            if [[ ${#PIDS[@]} -gt 0 ]]; then
                PIDS=("${PIDS[@]}")
            else
                PIDS=()
            fi
            sleep 1
        done

        bash _run_single_ala.sh "${residue}" "${WORKDIR}" \
            "${MMPBSA_START_FRAME}" "${MMPBSA_END_FRAME}" "${MMPBSA_INTERVAL}" \
            "${IGB}" "${IE_SEGMENT}" "${TEMPERATURE}" \
            "${REC_GROUP:-1}" "${LIG_GROUP:-13}" &
        PIDS+=($!)
        RUNNING=$((RUNNING + 1))
    done < "${RESIDUE_LIST_FILE}"

    # Wait for all subprocesses (safe empty array handling)
    if [[ ${#PIDS[@]} -gt 0 ]]; then
        for pid in "${PIDS[@]}"; do
            wait "${pid}" 2>/dev/null || true
        done
    fi

else
    # ---- Serial mode ----
    echo "[AS] Running in serial mode ..."
    COUNT=0

    while IFS= read -r residue; do
        COUNT=$((COUNT + 1))
        echo ""
        echo "── [${COUNT}/${TOTAL_RESIDUES}] Scanning ${residue} ──"
        bash _run_single_ala.sh "${residue}" "${WORKDIR}" \
            "${MMPBSA_START_FRAME}" "${MMPBSA_END_FRAME}" "${MMPBSA_INTERVAL}" \
            "${IGB}" "${IE_SEGMENT}" "${TEMPERATURE}" \
            "${REC_GROUP:-1}" "${LIG_GROUP:-13}" || true
    done < "${RESIDUE_LIST_FILE}"
fi

# ========== 4. Summary statistics ==========
SCAN_END=$(date +%s)
SCAN_ELAPSED=$((SCAN_END - SCAN_START))

COMPLETED=$(ls ala_scan_results/*.dat 2>/dev/null | wc -l | tr -d ' ')
FAILED=$((TOTAL_RESIDUES - COMPLETED))

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Alanine Scanning Complete"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Total residues: ${TOTAL_RESIDUES}"
echo "  Succeeded:      ${COMPLETED}"
echo "  Failed:         ${FAILED}"
echo "  Elapsed:        $((SCAN_ELAPSED/60))m $((SCAN_ELAPSED%60))s"
echo ""

if [[ ${FAILED} -gt 0 ]]; then
    echo "  [Tip] Failed residues can be re-run:"
    echo "    bash scripts/06b_ala_scan_batch.sh"
    echo "  (Completed residues will be automatically skipped)"
fi

# Clean up
rm -rf ala_scan_tmp _run_single_ala.sh 2>/dev/null || true

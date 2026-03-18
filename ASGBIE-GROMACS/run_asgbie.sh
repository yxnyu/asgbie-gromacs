#!/usr/bin/env bash
# ============================================================
#   ASGBIE-GROMACS  One-Click Run Script
#   ================================================
#   Fully replaces commercial AMBER ASGBIE workflow
#
#   Technical stack:
#     MD engine:       GROMACS (free & open-source)
#     Free energy:     gmx_MMPBSA + AmberTools (free)
#     ASGBIE:          Alanine Scanning + GB + Interaction Entropy
#
#   Usage:
#     1. Edit config.sh to set your system parameters
#     2. Place complex PDB in this directory
#     3. conda activate autodock_vina
#     4. bash run_asgbie.sh [--skip-md] [--only-asgbie] [--step N]
# ============================================================
set -euo pipefail

PROJ_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PROJ_DIR}/config.sh"

# ---------- Argument Parsing ----------
SKIP_MD=false
ONLY_ASGBIE=false
START_STEP=1

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-md)     SKIP_MD=true; START_STEP=6; shift ;;
        --only-asgbie) ONLY_ASGBIE=true; START_STEP=6; shift ;;
        --step)        START_STEP="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: bash run_asgbie.sh [options]"
            echo ""
            echo "Options:"
            echo "  --skip-md       Skip MD, start directly from ASGBIE calculation"
            echo "  --only-asgbie   Only run ASGBIE calculation (requires existing prod.tpr/xtc)"
            echo "  --step N        Start from step N (1-8)"
            echo "  -h, --help      Show help"
            echo ""
            echo "Steps:"
            echo "  1. System preparation (pdb2gmx)"
            echo "  2. Solvation + ions"
            echo "  3. Energy minimization"
            echo "  4. NVT + NPT equilibration"
            echo "  5. Production MD"
            echo "  6. ASGBIE calculation (core)"
            echo "  7. Result collection"
            echo "  8. Analysis & plotting"
            exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1"; exit 1 ;;
    esac
done

# ---------- Environment Check ----------
echo "╔══════════════════════════════════════════╗"
echo "║     ASGBIE-GROMACS  One-Click Workflow   ║"
echo "║     Replacing AMBER · 100% Free          ║"
echo "╚══════════════════════════════════════════╝"
echo ""
echo "Configuration:"
echo "  Input PDB:    ${INPUT_PDB}"
echo "  Force field:  ${FORCEFIELD}"
echo "  Temperature:  ${TEMPERATURE} K"
echo "  MD length:    $(echo "scale=1; ${PROD_STEPS}*0.002/1000" | bc) ns"
echo "  GB model:     igb=${IGB}"
echo "  Work dir:     ${WORKDIR}"
echo "  Start step:   Step ${START_STEP}"
echo ""

# Check essential tools
for cmd in gmx gmx_MMPBSA python3; do
    if ! command -v ${cmd} &>/dev/null; then
        echo "[ERROR] ${cmd} not found. Please activate the conda environment:"
        echo "  conda activate autodock_vina"
        exit 1
    fi
done
echo "[OK] GROMACS $(gmx --version 2>&1 | head -1 | awk '{print $NF}')"
echo "[OK] gmx_MMPBSA $(gmx_MMPBSA --version 2>&1 | head -1 || echo 'installed')"
echo ""

# Check input file
if [[ ${START_STEP} -le 1 ]] && [[ ! -f "${PROJ_DIR}/${INPUT_PDB}" ]]; then
    echo "[ERROR] Input file not found: ${PROJ_DIR}/${INPUT_PDB}"
    echo "  Please place the protein-ligand complex PDB in the project directory"
    exit 1
fi

# ---------- Timer ----------
TOTAL_START=$(date +%s)

run_step() {
    local step_num=$1
    local script_name=$2
    local desc=$3

    if [[ ${step_num} -lt ${START_STEP} ]]; then
        return
    fi

    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Step ${step_num}/8: ${desc}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    local step_start=$(date +%s)
    bash "${PROJ_DIR}/scripts/${script_name}"
    local step_end=$(date +%s)
    local elapsed=$((step_end - step_start))
    echo "  Elapsed: $((elapsed/60))m $((elapsed%60))s"
}

# ---------- Execute Workflow ----------
if [[ "${ONLY_ASGBIE}" == "false" ]]; then
    run_step 1 "01_prepare.sh"    "System Preparation"
    run_step 2 "02_solvate.sh"    "Solvation + Ions"
    run_step 3 "03_minimize.sh"   "Energy Minimization"
    run_step 4 "04_equilibrate.sh" "NVT + NPT Equilibration"
    run_step 5 "05_production.sh" "Production MD"
fi

run_step 6 "06_asgbie.sh"   "ASGBIE Free Energy Calculation"
run_step 7 "07_collect.sh"  "Result Collection"

# Step 8: Python plotting
if [[ ${START_STEP} -le 8 ]]; then
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Step 8/8: Analysis & Plotting"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    cd "${PROJ_DIR}"
    python3 scripts/08_analysis.py
fi

# ---------- Summary ----------
TOTAL_END=$(date +%s)
TOTAL_ELAPSED=$((TOTAL_END - TOTAL_START))

echo ""
echo "╔══════════════════════════════════════════╗"
echo "║       ASGBIE-GROMACS Complete             ║"
echo "╚══════════════════════════════════════════╝"
echo ""
echo "  Total time: $((TOTAL_ELAPSED/3600))h $((TOTAL_ELAPSED%3600/60))m $((TOTAL_ELAPSED%60))s"
echo ""
echo "  Results directory: ${RESULTS_DIR}/"
echo "  Key files:"
echo "    - FINAL_RESULTS_ASGBIE.dat  (ASGBIE binding free energy)"
echo "    - FINAL_RESULTS_ASGBIE.csv  (Detailed energy decomposition)"
echo "    - SUMMARY.txt               (Analysis summary)"
echo "    - plot_*.png                 (Plots)"
echo ""
echo "  Tip: Use gmx_MMPBSA_ana for interactive result viewing"
echo "    gmx_MMPBSA_ana -p ${WORKDIR}/asgbie_calc/_GMXMMPBSA_info"

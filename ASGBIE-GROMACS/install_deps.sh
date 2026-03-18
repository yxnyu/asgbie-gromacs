#!/usr/bin/env bash
# ============================================================
# ASGBIE-GROMACS Dependency Installation / Verification Script
# Uses the existing autodock_vina environment by default
# ============================================================
set -euo pipefail

ENV_NAME="${1:-autodock_vina}"

echo "============================================"
echo "  ASGBIE-GROMACS Environment Verification"
echo "  Environment: ${ENV_NAME}"
echo "============================================"

# ---------- Verify existing environment ----------
echo ""
echo "[Check] Core dependencies:"

check_pkg() {
    local name=$1
    if conda list -n ${ENV_NAME} 2>/dev/null | grep -qi "^${name} " ; then
        local ver=$(conda list -n ${ENV_NAME} 2>/dev/null | grep -i "^${name} " | awk '{print $2}')
        echo "  ✓ ${name} ${ver}"
        return 0
    else
        echo "  ✗ ${name} — not installed"
        return 1
    fi
}

ALL_OK=true
for pkg in gromacs ambertools python numpy pandas matplotlib scipy seaborn; do
    check_pkg ${pkg} || ALL_OK=false
done

# Check gmx_MMPBSA separately (pip package)
if conda run -n ${ENV_NAME} pip show gmx_MMPBSA &>/dev/null 2>&1; then
    ver=$(conda run -n ${ENV_NAME} pip show gmx_MMPBSA 2>/dev/null | grep Version | awk '{print $2}')
    echo "  ✓ gmx_MMPBSA ${ver}"
else
    echo "  ✗ gmx_MMPBSA — not installed"
    ALL_OK=false
fi

echo ""

if [ "${ALL_OK}" = true ]; then
    echo "[OK] All dependencies are ready, no additional installation needed"
    echo ""
    echo "Usage:"
    echo "  conda activate ${ENV_NAME}"
    echo "  bash run_asgbie.sh"
else
    echo "[WARN] Some dependencies are missing, attempting installation ..."
    INSTALLER=$(command -v mamba || command -v conda)
    ${INSTALLER} install -n ${ENV_NAME} -y -c conda-forge gromacs ambertools numpy pandas matplotlib scipy seaborn
    conda run -n ${ENV_NAME} pip install gmx_MMPBSA
    echo "[OK] Installation complete"
fi

#!/usr/bin/env bash
# ============================================================
# Step 3: Energy Minimization
# Replaces AMBER's min1/min2 steps
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "===== [Step 3] Energy Minimization ====="
cd "${WORKDIR}"

# ---------- Prepare MDP ----------
cp "${SCRIPT_DIR}/../mdp/em.mdp" em.mdp
sed -i.bak "s/__EM_STEPS__/${EM_STEPS}/g" em.mdp

# ---------- Generate tpr ----------
echo "[3.1] Generating em.tpr ..."
gmx grompp -f em.mdp \
    -c complex_ions.gro \
    -p topol.top \
    -o em.tpr \
    -maxwarn 5

# ---------- Run minimization ----------
echo "[3.2] Running energy minimization ..."
gmx mdrun -v -deffnm em ${GMX_GPU} \
    -ntmpi ${GMX_NTMPI} -ntomp ${GMX_NTOMP}

# ---------- Check convergence ----------
echo "[3.3] Checking potential energy convergence ..."
echo "Potential" | gmx energy -f em.edr -o em_potential.xvg

echo "[Step 3] Done ✓"
echo "  Final potential energy: $(tail -1 em_potential.xvg | awk '{print $2}') kJ/mol"

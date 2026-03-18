#!/usr/bin/env bash
# ============================================================
# Step 4: NVT + NPT Equilibration
# Replaces AMBER's heat1/heat2/back/Calpha steps
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "===== [Step 4] Equilibration ====="
cd "${WORKDIR}"

# ========== NVT ==========
echo "[4.1] NVT equilibration (${TEMPERATURE} K) ..."
cp "${SCRIPT_DIR}/../mdp/nvt.mdp" nvt.mdp
sed -i.bak "s/__NVT_STEPS__/${NVT_STEPS}/g; s/__TEMPERATURE__/${TEMPERATURE}/g" nvt.mdp

gmx grompp -f nvt.mdp -c em.gro -r em.gro \
    -p topol.top -o nvt.tpr -maxwarn 5

gmx mdrun -v -deffnm nvt ${GMX_GPU} \
    -ntmpi ${GMX_NTMPI} -ntomp ${GMX_NTOMP}

echo "  NVT temperature check:"
echo "Temperature" | gmx energy -f nvt.edr -o nvt_temp.xvg
echo "  Average temperature: $(tail -5 nvt_temp.xvg | awk '{s+=$2;n++}END{printf "%.1f K",s/n}')"

# ========== NPT ==========
echo "[4.2] NPT equilibration (1 bar) ..."
cp "${SCRIPT_DIR}/../mdp/npt.mdp" npt.mdp
sed -i.bak "s/__NPT_STEPS__/${NPT_STEPS}/g; s/__TEMPERATURE__/${TEMPERATURE}/g" npt.mdp

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro \
    -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 5

gmx mdrun -v -deffnm npt ${GMX_GPU} \
    -ntmpi ${GMX_NTMPI} -ntomp ${GMX_NTOMP}

echo "  NPT density check:"
echo "Density" | gmx energy -f npt.edr -o npt_density.xvg
echo "  Average density: $(tail -5 npt_density.xvg | awk '{s+=$2;n++}END{printf "%.1f kg/m³",s/n}')"

echo "[Step 4] Done ✓"

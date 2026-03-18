#!/usr/bin/env bash
# ============================================================
# Step 5: Production MD
# Replaces AMBER's pmemd.cuda prod step
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "===== [Step 5] Production MD ====="
cd "${WORKDIR}"

# ---------- Prepare MDP ----------
cp "${SCRIPT_DIR}/../mdp/prod.mdp" prod.mdp
sed -i.bak \
    "s/__PROD_STEPS__/${PROD_STEPS}/g; \
     s/__PROD_NSTXOUT__/${PROD_NSTXOUT}/g; \
     s/__TEMPERATURE__/${TEMPERATURE}/g" prod.mdp

# ---------- Generate tpr ----------
echo "[5.1] Generating prod.tpr ..."
gmx grompp -f prod.mdp -c npt.gro -r npt.gro \
    -t npt.cpt -p topol.top -o prod.tpr -maxwarn 5

# ---------- Run production MD ----------
TOTAL_NS=$(echo "scale=1; ${PROD_STEPS} * 0.002 / 1000" | bc)
echo "[5.2] Running production MD (${TOTAL_NS} ns) ..."
echo "  Note: This step takes a long time, please be patient"

gmx mdrun -v -deffnm prod ${GMX_GPU} \
    -ntmpi ${GMX_NTMPI} -ntomp ${GMX_NTOMP}

# ---------- RMSD Analysis ----------
echo "[5.3] Computing RMSD ..."
# Backbone RMSD
echo -e "Backbone\nBackbone" | gmx rms -s prod.tpr -f prod.xtc \
    -o rmsd_backbone.xvg -tu ns

# Protein RMSD
echo -e "Protein\nProtein" | gmx rms -s prod.tpr -f prod.xtc \
    -o rmsd_protein.xvg -tu ns

echo "[Step 5] Done ✓"
echo "  Trajectory: prod.xtc ($(du -h prod.xtc | cut -f1))"
echo "  RMSD: rmsd_backbone.xvg, rmsd_protein.xvg"

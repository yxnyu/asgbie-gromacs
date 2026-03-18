#!/usr/bin/env bash
# ============================================================
# Step 2: Solvation + Ion Addition
# Replaces AMBER's tLEaP solvation step
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "===== [Step 2] Solvation & Ions ====="
cd "${WORKDIR}"

# ---------- Define box ----------
echo "[2.1] Defining simulation box (${BOX_TYPE}, margin ${BOX_DIST} nm) ..."
gmx editconf -f complex_processed.gro \
    -o complex_box.gro \
    -c -d ${BOX_DIST} \
    -bt ${BOX_TYPE}

# ---------- Add solvent ----------
echo "[2.2] Adding water molecules ..."
gmx solvate -cp complex_box.gro \
    -cs spc216.gro \
    -o complex_solv.gro \
    -p topol.top

# ---------- Add ions ----------
echo "[2.3] Adding ions (concentration ${ION_CONC} M) ..."

# Generate a simple tpr for genion
SCRIPT_DIR_ABS="$(cd "${SCRIPT_DIR}" && pwd)"
cp "${SCRIPT_DIR_ABS}/../mdp/em.mdp" em_tmp.mdp
sed -i.bak "s/__EM_STEPS__/${EM_STEPS}/g" em_tmp.mdp

gmx grompp -f em_tmp.mdp \
    -c complex_solv.gro \
    -p topol.top \
    -o ions.tpr \
    -maxwarn 5

# SOL group used to replace water molecules with ions
echo "SOL" | gmx genion -s ions.tpr \
    -o complex_ions.gro \
    -p topol.top \
    -pname ${PNAME} -nname ${NNAME} \
    -conc ${ION_CONC} -neutral

rm -f em_tmp.mdp em_tmp.mdp.bak ions.tpr

echo "[Step 2] Done ✓"
echo "  Output: complex_ions.gro, topol.top (updated)"

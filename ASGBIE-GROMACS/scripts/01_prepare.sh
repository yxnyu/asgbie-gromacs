#!/usr/bin/env bash
# ============================================================
# Step 1: System Preparation — PDB → GROMACS Topology
# Replaces AMBER's tLEaP
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "===== [Step 1] System Preparation ====="
mkdir -p "${WORKDIR}" && cd "${WORKDIR}"

# Copy input PDB
cp "${SCRIPT_DIR}/../${INPUT_PDB}" ./complex_raw.pdb

# ---------- Separate receptor and ligand chains ----------
echo "[1.1] Separating receptor (${RECEPTOR_CHAINS}) and ligand (${LIGAND_CHAINS}) chains ..."

python3 << 'PYEOF'
import sys, os
receptor_chains = os.environ.get("RECEPTOR_CHAINS", "A").split()
ligand_chains   = os.environ.get("LIGAND_CHAINS", "B").split()

rec_lines, lig_lines, com_lines = [], [], []
with open("complex_raw.pdb") as f:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            chain = line[21]
            com_lines.append(line)
            if chain in receptor_chains:
                rec_lines.append(line)
            elif chain in ligand_chains:
                lig_lines.append(line)

for fname, lines in [("receptor.pdb", rec_lines), ("ligand.pdb", lig_lines), ("complex.pdb", com_lines)]:
    with open(fname, "w") as out:
        out.writelines(lines)
        out.write("END\n")
    print(f"  Written {fname}: {len(lines)} lines")
PYEOF

# ---------- Generate GROMACS topology ----------
echo "[1.2] Running pdb2gmx (force field: ${FORCEFIELD}, water model: ${WATER_MODEL}) ..."
gmx pdb2gmx -f complex.pdb \
    -o complex_processed.gro \
    -p topol.top \
    -i posre.itp \
    -ff ${FORCEFIELD} \
    -water ${WATER_MODEL} \
    -ignh -merge all <<< "1"

echo "[1.3] Generating index file ..."
# Create receptor + ligand index groups
echo -e "chain ${RECEPTOR_CHAINS}\nchain ${LIGAND_CHAINS}\nq" | \
    gmx make_ndx -f complex_processed.gro -o index.ndx 2>/dev/null || true

echo "[Step 1] Done ✓"
echo "  Output: complex_processed.gro, topol.top, index.ndx"

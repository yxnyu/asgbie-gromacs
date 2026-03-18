#!/usr/bin/env bash
# ============================================================
# ASGBIE-GROMACS User Configuration
# Modify the variables below to fit your system
# ============================================================

# ---------- Environment ----------
CONDA_ENV="autodock_vina"

# ---------- Input Files ----------
# Protein-ligand complex PDB (ensure water and heteroatoms are removed)
INPUT_PDB="complex.pdb"

# ---------- Force Field & Water Model ----------
# Options: charmm27, amber99sb-ildn, amber03, oplsaa
FORCEFIELD="charmm27"
WATER_MODEL="tip3p"

# ---------- System Definition ----------
# Receptor chain IDs (space-separated)
RECEPTOR_CHAINS="A"
# Ligand chain IDs (space-separated)
LIGAND_CHAINS="B"

# ---------- Box Parameters ----------
BOX_TYPE="dodecahedron"     # cubic / dodecahedron / octahedron
BOX_DIST="1.2"              # Distance from protein to box edge (nm)

# ---------- Ions ----------
ION_CONC="0.15"             # NaCl concentration (mol/L)
PNAME="NA"                  # Cation name
NNAME="CL"                  # Anion name

# ---------- MD Parameters ----------
EM_STEPS="50000"            # Energy minimization max steps
NVT_STEPS="50000"           # NVT equilibration (100 ps @ dt=0.002)
NPT_STEPS="50000"           # NPT equilibration (100 ps @ dt=0.002)
PROD_STEPS="5000000"        # Production MD (10 ns @ dt=0.002)
PROD_NSTXOUT="5000"         # Trajectory output frequency (every 10 ps)
TEMPERATURE="300"           # Temperature (K)

# ---------- GPU Acceleration ----------
GMX_GPU=""                   # GPU: "-nb gpu", CPU only: ""
GMX_NTMPI="1"               # MPI ranks
GMX_NTOMP="8"               # OpenMP threads

# ---------- ASGBIE Calculation Parameters ----------
# Frame range for MM/GBSA analysis
MMPBSA_START_FRAME="1"
MMPBSA_END_FRAME="100"
MMPBSA_INTERVAL="1"

# GB model (1=HCT, 2=OBC1, 5=OBC2, recommended igb=5)
IGB="5"

# Interaction Entropy
IE_SEGMENT="25"             # IE segment count

# ---------- Alanine Scanning ----------
# Enable/disable alanine scanning (true/false)
ALA_SCAN_ENABLED="true"

# Residue list for mutation (comma-separated, empty = auto-detect interface residues)
# Format: "RESNAME:RESNUM" e.g., "ARG:25,GLU:30,TYR:45"
ALA_SCAN_RESIDUES=""

# Max residues to scan during auto-detection (prevents excessive runtime)
ALA_SCAN_MAX_RESIDUES="50"

# Parallel jobs (1=serial, >1=parallel, recommended: half of CPU cores)
# Requires GNU Parallel: conda install -c conda-forge parallel
ALA_SCAN_PARALLEL="1"

# ---------- Index Group Numbers ----------
# Receptor and ligand group numbers from gmx make_ndx
# Before first run, check: gmx make_ndx -f prod.tpr to see group numbers
REC_GROUP="1"                # Receptor group number
LIG_GROUP="13"               # Ligand group number

# ---------- Output ----------
WORKDIR="$(pwd)/asgbie_work"
RESULTS_DIR="$(pwd)/asgbie_results"

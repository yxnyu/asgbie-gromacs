# ASGBIE-GROMACS Technical Report

## Complete Migration from Commercial AMBER to a Fully Free Open-Source Solution

---

## 1. Background & Motivation

### 1.1 ASGBIE Method Overview

**ASGBIE** (Alanine Scanning Generalized Born with Interaction Entropy) is a protein-ligand binding free energy calculation method developed by the John Z.H. Zhang group, combining three techniques:

| Component | Full Name | Function |
|-----------|-----------|----------|
| **AS** | Alanine Scanning | Mutate interface residues to alanine one by one, evaluate per-residue contribution |
| **GB** | Generalized Born | Implicit solvent model for efficient solvation free energy calculation |
| **IE** | Interaction Entropy | Estimate entropy contribution from energy fluctuations, replacing traditional Normal Mode or QH methods |

### 1.2 Why Migrate to GROMACS

The original ASGBIE workflow is entirely based on the **AMBER** software suite (pmemd.cuda + MMPBSA.py), with the following issues:

| Issue | Description |
|-------|-------------|
| **Commercial license** | Full AMBER requires payment (academic ~$500, commercial even more) |
| **pmemd.cuda** | GPU-accelerated MD engine only available in paid version |
| **Closed ecosystem** | Scripts highly dependent on AMBER-specific formats (prmtop, inpcrd) |

### 1.3 Solution

This project uses a **100% free open-source toolchain** as replacement:

| AMBER Component | Replacement | License |
|-----------------|-------------|---------|
| pmemd.cuda (MD) | **GROMACS** 2025.4 | LGPL |
| tLEaP (modeling) | **gmx pdb2gmx** | LGPL |
| ff15ipq + Lipid21 | **CHARMM36m** (or AMBER ported) | Free |
| MMPBSA.py | **gmx_MMPBSA** + AmberTools | GPL / Free |
| Alanine scanning | gmx_MMPBSA `&alanine_scanning` | Free |
| Interaction entropy | gmx_MMPBSA `interaction_entropy=1` | Free |

**Key technical point**: `gmx_MMPBSA` internally converts GROMACS topology (.top) and trajectory (.xtc) to AMBER format, then calls the **AmberTools** (free open-source version) MMPBSA.py for computation. This completely bypasses AMBER's commercial license restriction.

---

## 2. Technical Architecture

### 2.1 Workflow Comparison

```
Original AMBER workflow (requires commercial license):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
tLEaP → PACKMOL-Memgen → pmemd.cuda (min/heat/prod) → MMPBSA.py → Results
  │         │                │                            │
  └─ $$$  ──┘────── $$$$ ───┘──────────── $$$ ───────────┘

GROMACS workflow (all free):
━━━━━━━━━━━━━━━━━━━━━━━
pdb2gmx → solvate → gmx mdrun (GPU) → gmx_MMPBSA → Results
  │          │            │                │
  └─ FREE ──┘──── FREE ──┘──── FREE ─────┘
```

### 2.2 8-Step Workflow

| Step | Script | Function | AMBER Equivalent |
|------|--------|----------|-----------------|
| 1 | `01_prepare.sh` | PDB → GRO/TOP | tLEaP |
| 2 | `02_solvate.sh` | Solvation + ions | tLEaP addions |
| 3 | `03_minimize.sh` | Energy minimization | pmemd min1/min2 |
| 4 | `04_equilibrate.sh` | NVT + NPT equilibration | pmemd heat/back/Calpha |
| 5 | `05_production.sh` | Production MD + RMSD | pmemd prod + cpptraj |
| 6 | `06_asgbie.sh` | **ASGBIE core calculation** | MMPBSA.py |
| 6b | `06b_ala_scan_batch.sh` | **Batch alanine scanning** | Manual loop scripts |
| 7 | `07_collect.sh` | Result collection | Custom scripts |
| 8 | `08_analysis.py` | Plot generation | RMSD.py |

---

## 3. Method Details

### 3.1 MM/GBSA Binding Free Energy

Binding free energy decomposition:

```
ΔG_bind = ΔH - TΔS
        = ΔE_MM + ΔG_solv - TΔS

Where:
  ΔE_MM  = ΔE_vdw + ΔE_ele          (Molecular mechanics energy)
  ΔG_solv = ΔG_GB + ΔG_SA            (Solvation free energy)
  TΔS    = Interaction Entropy (IE)   (Entropy contribution)
```

### 3.2 Generalized Born (GB) Model

This solution uses **igb=5** (OBC2 model):

- High computational efficiency (10-100x faster than PB)
- Good accuracy for protein-protein/peptide interactions
- Salt concentration simulated via `saltcon=0.15` for physiological conditions

### 3.3 Interaction Entropy (IE)

Traditional MM/GBSA uses Normal Mode Analysis (NMA) or Quasi-Harmonic (QH) for entropy:
- NMA: Computationally extremely expensive
- QH: Requires extensive conformational sampling

**IE method** based on energy fluctuation theory:

```
-TΔS_IE = kT × ln⟨exp(β × ΔE_interaction)⟩

Where β = 1/kT, ΔE_interaction = E_interaction - ⟨E_interaction⟩
```

Advantages:
- No additional computation needed, extracted directly from MD trajectory
- Computational overhead is nearly zero
- Accuracy superior to traditional NMA/QH

### 3.4 Alanine Scanning

For each non-alanine/glycine residue at the interface, perform in silico mutation:

```
ΔΔG_bind = ΔG_bind(mutant) - ΔG_bind(wildtype)
```

- ΔΔG > 0: The residue positively contributes to binding ("hotspot" residue)
- ΔΔG < 0: The residue is unfavorable for binding

### 3.5 Batch Alanine Scanning Improvement

The original `gmx_MMPBSA` `&alanine_scanning` only supports single-residue mutation per run. This project implements a batch approach:

- **Per-residue loop**: Each residue is scanned in an independent gmx_MMPBSA run
- **Resume support**: Completed residues are automatically skipped on re-run
- **Parallel execution**: GNU Parallel or background processes for multi-core acceleration
- **Fault tolerance**: One residue failure does not affect others

---

## 4. Software Dependencies

All dependencies can be obtained for free via conda:

| Software | Version | Purpose | Installation |
|----------|---------|---------|-------------|
| GROMACS | 2025.4 | MD simulation engine | `conda install gromacs` |
| AmberTools | 24.8 | MMPBSA backend | `conda install ambertools` |
| gmx_MMPBSA | 1.6.4 | GROMACS→AMBER bridge | `pip install gmx_MMPBSA` |
| Python | 3.11 | Scripts/analysis | conda |
| NumPy | 2.4 | Numerical computation | conda |
| Pandas | 3.0 | Data processing | conda |
| Matplotlib | 3.10 | Plotting | conda |
| Seaborn | 0.13 | Statistical plotting | conda |

**Total license cost: $0**

---

## 5. Usage Guide

### 5.1 Quick Start

```bash
# 1. Activate environment
conda activate autodock_vina

# 2. Enter project directory
cd ASGBIE-GROMACS

# 3. Edit configuration (set PDB name, chain IDs, etc.)
vim config.sh

# 4. Place your complex PDB
cp /path/to/your/complex.pdb .

# 5. One-click run
bash run_asgbie.sh
```

### 5.2 Advanced Usage

```bash
# Only run ASGBIE calculation (MD trajectory already exists)
bash run_asgbie.sh --only-asgbie

# Start from a specific step
bash run_asgbie.sh --step 4    # Start from NPT equilibration

# Skip MD, directly calculate free energy
bash run_asgbie.sh --skip-md
```

### 5.3 Key Configuration Items

In `config.sh`:

```bash
# Your protein-ligand complex
INPUT_PDB="complex.pdb"
RECEPTOR_CHAINS="A"       # Receptor chains
LIGAND_CHAINS="B"         # Ligand chains

# MD parameters
PROD_STEPS="5000000"      # 10 ns (can increase to 50M = 100 ns)
TEMPERATURE="300"         # K

# ASGBIE parameters
IGB="5"                   # GB model
IE_SEGMENT="25"           # IE segment count

# Alanine scanning
ALA_SCAN_PARALLEL="4"     # Parallel jobs
ALA_SCAN_RESIDUES=""      # Empty = auto-detect, or specify manually
```

---

## 6. gmx_MMPBSA Input File Explained

Core input file `mmpbsa_asgbie.in`:

```
&general
  interaction_entropy = 1,    # ← Enable IE (the "IE" part of ASGBIE)
  ie_segment = 25,            # IE segment count
  temperature = 300,          # Temperature
/

&gb
  igb = 5,                    # ← OBC2 GB model (the "GB" part of ASGBIE)
  saltcon = 0.15,             # Ionic strength
/

&alanine_scanning
  mutant = "ALA",             # ← Alanine scanning (the "AS" part of ASGBIE)
/
```

The combination of these three namelists = **complete ASGBIE method**.

---

## 7. Differences from Original AMBER Version

| Aspect | AMBER Original | GROMACS Version |
|--------|---------------|-----------------|
| MD engine | pmemd.cuda | gmx mdrun (GPU) |
| Force field | ff15ipq + Lipid21 | CHARMM36m (or AMBER ported) |
| Water model | SPC/Eb | TIP3P (or SPC/E) |
| MM/GBSA | MMPBSA.py (direct) | gmx_MMPBSA (bridge) |
| Alanine scanning | Manual loop scripts | gmx_MMPBSA built-in + batch loop |
| Interaction entropy | Custom calculation | gmx_MMPBSA built-in |
| Membrane protein | PACKMOL-Memgen | CHARMM-GUI (extensible) |
| Cost | ~$500+ | **$0** |
| GPU acceleration | Paid version only | **Free** |

### 7.1 Accuracy Expectations

- MM/GBSA: GROMACS and AMBER use the same MMPBSA.py backend, accuracy is identical
- Force field difference: CHARMM36m vs ff15ipq may introduce ~0.5-1.0 kcal/mol systematic bias
- IE method: Same implementation, results should be identical
- Alanine scanning: Same calculation method, residue ranking (relative values) highly consistent

---

## 8. References

1. Zhou M, et al. "An efficient approach to the accurate prediction of mutational effects in antigen binding to the MHC1." *Molecules* (2024). — **ASGBIE core paper**

2. Valdés-Tresanco MS, et al. "gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS." *J. Chem. Theory Comput.* (2021). — **gmx_MMPBSA method paper**

3. Zhang J, et al. "Combined antibodies evusheld against SARS-CoV-2 omicron variants." *J. Chem. Inf. Model.* (2023). — **ASGBIE application**

4. Li J, et al. "Significantly enhancing human antibody affinity via deep learning and ASGBIE." *Briefings in Bioinformatics* (2025). — **ASGBIE + AI**

5. Abraham MJ, et al. "GROMACS: High performance molecular simulations." *SoftwareX* (2015). — **GROMACS**

6. GitHub: [xxttzhang/MD-ASGBIE](https://github.com/xxttzhang/MD-ASGBIE) — **Original AMBER scripts**

---

## 9. Project File Structure

```
ASGBIE-GROMACS/
├── run_asgbie.sh              # One-click main script
├── config.sh                  # User configuration (edit this)
├── install_deps.sh            # Environment verification/installation
├── complex.pdb                # ← Your input file goes here
│
├── scripts/
│   ├── 01_prepare.sh          # PDB → GROMACS topology
│   ├── 02_solvate.sh          # Solvation + ions
│   ├── 03_minimize.sh         # Energy minimization
│   ├── 04_equilibrate.sh      # NVT + NPT equilibration
│   ├── 05_production.sh       # Production MD
│   ├── 06_asgbie.sh           # ★ ASGBIE core calculation
│   ├── 06b_ala_scan_batch.sh  # Batch alanine scanning (per-residue loop)
│   ├── 07_collect.sh          # Result collection
│   └── 08_analysis.py         # Analysis & plotting
│
├── mdp/
│   ├── em.mdp                 # Energy minimization parameters
│   ├── nvt.mdp                # NVT equilibration parameters
│   ├── npt.mdp                # NPT equilibration parameters
│   ├── prod.mdp               # Production MD parameters
│   └── mmpbsa_asgbie.in       # ASGBIE input template
│
├── docs/
│   └── REPORT.md              # This report
│
└── example/                   # Example files
```

---

*Report generated: 2026-03-18*
*ASGBIE-GROMACS v1.0*

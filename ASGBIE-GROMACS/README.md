# ASGBIE-GROMACS

**100% Free & Open-Source Protein-Ligand Binding Free Energy Workflow**

Fully replaces the commercial AMBER suite using GROMACS + gmx_MMPBSA + AmberTools to implement the ASGBIE (Alanine Scanning + Generalized Born + Interaction Entropy) method.

---

## Features

- **Zero cost**: All tools are free and open-source (GROMACS + AmberTools + gmx_MMPBSA)
- **GPU accelerated**: GROMACS natively supports GPU-accelerated MD simulations
- **One-click run**: `bash run_asgbie.sh` from PDB to final results
- **Batch alanine scanning**: Per-residue loop calling gmx_MMPBSA with resume support and parallel acceleration
- **Automated analysis**: Generates RMSD, energy decomposition, hotspot residue waterfall/heatmap plots

## ASGBIE Method

| Component | Full Name | Function |
|-----------|-----------|----------|
| **AS** | Alanine Scanning | Mutate interface residues to alanine one by one, evaluate per-residue contribution |
| **GB** | Generalized Born | Implicit solvent MM/GBSA binding free energy (igb=5, OBC2) |
| **IE** | Interaction Entropy | Estimate entropy contribution from energy fluctuations, replacing NMA/QH |

```
ΔG_bind = (ΔE_vdw + ΔE_ele + ΔG_GB + ΔG_SA) - TΔS_IE
```

## Quick Start

### 1. Environment Setup

```bash
# Create/verify conda environment
conda create -n asgbie python=3.11 -y
conda activate asgbie
conda install -c conda-forge gromacs ambertools numpy pandas matplotlib seaborn -y
pip install gmx_MMPBSA

# Or verify existing environment
bash install_deps.sh
```

### 2. Configure Your System

Edit `config.sh`:

```bash
INPUT_PDB="complex.pdb"    # Protein-ligand complex
RECEPTOR_CHAINS="A"         # Receptor chain ID
LIGAND_CHAINS="B"           # Ligand chain ID
PROD_STEPS="5000000"        # MD length (10 ns)
```

### 3. Run

```bash
# Place your complex PDB in the project root directory
cp /path/to/your/complex.pdb .

# One-click full workflow
bash run_asgbie.sh
```

### 4. View Results

```
asgbie_results/
├── FINAL_RESULTS_GBIE.dat     # GB+IE binding free energy
├── ala_scan_summary.csv        # Alanine scanning hotspot ranking
├── plot_rmsd.png               # RMSD plot
├── plot_ala_scan_waterfall.png  # Hotspot residue waterfall plot
├── plot_ala_scan_heatmap.png   # Hotspot residue heatmap
└── SUMMARY.txt                 # Analysis report
```

## 8-Step Workflow

| Step | Script | Function | AMBER Equivalent |
|------|--------|----------|-----------------|
| 1 | `01_prepare.sh` | PDB → GROMACS topology | tLEaP |
| 2 | `02_solvate.sh` | Solvation + ions | tLEaP addions |
| 3 | `03_minimize.sh` | Energy minimization | pmemd min |
| 4 | `04_equilibrate.sh` | NVT + NPT equilibration | pmemd heat |
| 5 | `05_production.sh` | Production MD + RMSD | pmemd prod |
| 6 | `06_asgbie.sh` | **ASGBIE core calculation** | MMPBSA.py |
| 6b | `06b_ala_scan_batch.sh` | **Batch alanine scanning** | Manual loop |
| 7 | `07_collect.sh` | Result collection | — |
| 8 | `08_analysis.py` | Plot generation | — |

## Advanced Usage

### Skip MD / Start from a Specific Step

```bash
# Already have MD trajectory, only run ASGBIE calculation
bash run_asgbie.sh --only-asgbie

# Start from NPT equilibration (previous steps already done)
bash run_asgbie.sh --step 4

# Skip MD, directly calculate free energy
bash run_asgbie.sh --skip-md
```

### Alanine Scanning Configuration

In `config.sh`:

```bash
# Enable/disable
ALA_SCAN_ENABLED="true"

# Manually specify residues (empty = auto-detect interface residues)
ALA_SCAN_RESIDUES="ARG:25,GLU:30,TYR:45,PHE:82"

# Parallel acceleration (requires GNU Parallel)
ALA_SCAN_PARALLEL="4"

# Maximum residues when auto-detecting
ALA_SCAN_MAX_RESIDUES="50"
```

### Improved Alanine Scanning Design

The original `gmx_MMPBSA` `&alanine_scanning` mutates all interface residues in a single run, which has issues:

| Problem | Solution in This Project |
|---------|------------------------|
| One residue failure aborts the entire run | **Per-residue independent execution**, failure does not affect others |
| Cannot specify a subset of residues | Supports `ALA_SCAN_RESIDUES` for precise specification |
| Cannot parallelize | Supports GNU Parallel multi-core parallelism |
| Cannot resume from interruption | Completed residues are automatically skipped |
| Scattered results hard to summarize | Automatically generates `ala_scan_summary.csv` ranking table |

```bash
# Re-run alanine scanning independently (completed residues auto-skipped)
bash scripts/06b_ala_scan_batch.sh

# Install GNU Parallel to enable parallelism
conda install -c conda-forge parallel
```

### Interactive Result Viewer

```bash
# gmx_MMPBSA built-in analysis tool
gmx_MMPBSA_ana -p asgbie_work/asgbie_calc/_GMXMMPBSA_info
```

## Project Structure

```
ASGBIE-GROMACS/
├── run_asgbie.sh              # One-click main script
├── config.sh                  # User configuration
├── install_deps.sh            # Environment verification/installation
├── complex.pdb                # ← Your input file goes here
│
├── scripts/
│   ├── 01_prepare.sh          # System preparation
│   ├── 02_solvate.sh          # Solvation
│   ├── 03_minimize.sh         # Energy minimization
│   ├── 04_equilibrate.sh      # Equilibration
│   ├── 05_production.sh       # Production MD
│   ├── 06_asgbie.sh           # ASGBIE core (GB+IE + alanine scanning)
│   ├── 06b_ala_scan_batch.sh  # Batch alanine scanning (per-residue loop)
│   ├── 07_collect.sh          # Result collection
│   └── 08_analysis.py         # Analysis & plotting
│
├── mdp/                       # GROMACS parameter files
│   ├── em.mdp                 # Energy minimization
│   ├── nvt.mdp                # NVT equilibration
│   ├── npt.mdp                # NPT equilibration
│   ├── prod.mdp               # Production MD
│   └── mmpbsa_asgbie.in       # MMPBSA input template
│
├── docs/
│   └── REPORT.md              # Technical report
│
└── test_hiv/                  # HIV protease test case
```

## Key Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `FORCEFIELD` | charmm27 | Force field (charmm27/amber99sb-ildn/oplsaa) |
| `PROD_STEPS` | 5000000 | MD steps (10 ns @ dt=0.002) |
| `IGB` | 5 | GB model (5=OBC2, recommended) |
| `IE_SEGMENT` | 25 | Interaction entropy segment count |
| `ALA_SCAN_PARALLEL` | 1 | Alanine scanning parallel jobs |
| `REC_GROUP` | 1 | Receptor index group number |
| `LIG_GROUP` | 13 | Ligand index group number |
| `BOX_TYPE` | dodecahedron | Simulation box type |
| `ION_CONC` | 0.15 | NaCl concentration (M) |

> **Note**: `REC_GROUP` and `LIG_GROUP` need to be adjusted for your system. Before the first run, execute `gmx make_ndx -f prod.tpr` to check group numbers.

## Software Dependencies

| Software | Version | Purpose | Installation |
|----------|---------|---------|-------------|
| GROMACS | ≥2022 | MD simulation | `conda install -c conda-forge gromacs` |
| AmberTools | ≥22 | MMPBSA backend | `conda install -c conda-forge ambertools` |
| gmx_MMPBSA | ≥1.6 | GROMACS↔AMBER bridge | `pip install gmx_MMPBSA` |
| Python | ≥3.9 | Scripts/analysis | conda |
| GNU Parallel | — | Parallel alanine scanning (optional) | `conda install -c conda-forge parallel` |

**Total license cost: $0**

## Comparison with Original AMBER Version

| Aspect | AMBER Original | GROMACS Version |
|--------|---------------|-----------------|
| MD engine | pmemd.cuda (paid) | gmx mdrun (free GPU) |
| Force field | ff15ipq | CHARMM36m / AMBER ported |
| MM/GBSA | MMPBSA.py | gmx_MMPBSA (same backend) |
| Alanine scanning | Manual loop | **Automated batch + parallel** |
| Resume support | None | **Completed residues auto-skipped** |
| Cost | ~$500+ | **$0** |

## References

1. Zhou M, et al. "An efficient approach to the accurate prediction of mutational effects in antigen binding to the MHC1." *Molecules* (2024). — ASGBIE method
2. Valdés-Tresanco MS, et al. "gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS." *JCTC* (2021). — gmx_MMPBSA
3. Li J, et al. "Significantly enhancing human antibody affinity via deep learning and ASGBIE." *Briefings in Bioinformatics* (2025). — ASGBIE + AI
4. Abraham MJ, et al. "GROMACS: High performance molecular simulations." *SoftwareX* (2015). — GROMACS
5. Original AMBER scripts: [xxttzhang/MD-ASGBIE](https://github.com/xxttzhang/MD-ASGBIE)

## License

MIT

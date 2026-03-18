<p align="center">
  <h1 align="center">ASGBIE-GROMACS</h1>
  <p align="center">
    <strong>100% Free & Open-Source Protein Binding Free Energy Workflow</strong>
  </p>
  <p align="center">
    <a href="#quick-start">Quick Start</a> •
    <a href="#8-step-workflow">Workflow</a> •
    <a href="#batch-alanine-scanning">Alanine Scanning</a> •
    <a href="docs/REPORT.md">Technical Report</a>
  </p>
  <p align="center">
    <img src="https://img.shields.io/badge/GROMACS-≥2022-blue?style=flat-square" alt="GROMACS">
    <img src="https://img.shields.io/badge/gmx__MMPBSA-≥1.6-green?style=flat-square" alt="gmx_MMPBSA">
    <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square" alt="License">
    <img src="https://img.shields.io/badge/Cost-$0-brightgreen?style=flat-square" alt="Cost">
  </p>
</p>

---

Fully replaces the **commercial AMBER** suite using **GROMACS + gmx_MMPBSA + AmberTools** to implement the ASGBIE method for protein–protein / protein–peptide binding free energy calculations.

## Why ASGBIE-GROMACS?

| | AMBER (Original) | **ASGBIE-GROMACS** |
|---|---|---|
| MD Engine | pmemd.cuda (**paid ~$500+**) | GROMACS mdrun (**free**, GPU) |
| MM/GBSA | MMPBSA.py | gmx_MMPBSA (same AmberTools backend) |
| Alanine Scanning | Manual scripts, single-run | **Batch per-residue loop + parallel** |
| Resume on Failure | ✗ | **✓ Auto-skip completed residues** |
| Total Cost | ~$500+ | **$0** |

## ASGBIE Method

**A**lanine **S**canning + **G**eneralized **B**orn + **I**nteraction **E**ntropy

```
ΔG_bind = (ΔE_vdw + ΔE_ele + ΔG_GB + ΔG_SA) − TΔS_IE
```

| Component | What it does |
|-----------|-------------|
| **AS** — Alanine Scanning | Mutate interface residues to Ala one-by-one → identify hotspot residues (ΔΔG > 1.0 kcal/mol) |
| **GB** — Generalized Born | Implicit solvent MM/GBSA binding energy (igb=5, OBC2 model) |
| **IE** — Interaction Entropy | Entropy from energy fluctuations — near-zero cost, replaces expensive NMA/QH |

## Features

- **One-click workflow** — `bash run_asgbie.sh` takes you from PDB to publication-ready results
- **Batch alanine scanning** — Per-residue independent runs with GNU Parallel support
- **Fault tolerant** — One residue failure doesn't kill the entire scan; re-run auto-skips completed residues
- **GPU accelerated** — GROMACS native GPU support, no license needed
- **Automated plots** — RMSD, energy decomposition, hotspot waterfall chart, heatmap

---

## Quick Start

### 1. Install Dependencies

```bash
conda create -n asgbie python=3.11 -y
conda activate asgbie
conda install -c conda-forge gromacs ambertools numpy pandas matplotlib seaborn -y
pip install gmx_MMPBSA

# Optional: parallel alanine scanning
conda install -c conda-forge parallel -y

# Or verify an existing environment
bash install_deps.sh
```

### 2. Configure

Edit `config.sh` for your system:

```bash
INPUT_PDB="complex.pdb"       # Your protein complex
RECEPTOR_CHAINS="A"            # Receptor chain(s)
LIGAND_CHAINS="B"              # Ligand/peptide chain(s)
PROD_STEPS="5000000"           # 10 ns production MD
TEMPERATURE="300"              # Kelvin
```

### 3. Run

```bash
cp /path/to/your/complex.pdb .
bash run_asgbie.sh
```

### 4. Results

```
asgbie_results/
├── FINAL_RESULTS_GBIE.dat        # GB + IE binding free energy
├── FINAL_RESULTS_GBIE.csv        # Detailed per-frame decomposition
├── ala_scan_summary.csv           # Hotspot residue ranking (ΔΔG)
├── ala_scan_details/              # Per-residue raw results
├── plot_rmsd.png                  # RMSD over time
├── plot_em_energy.png             # Energy minimization convergence
├── plot_equilibration.png         # NVT temperature + NPT density
├── plot_ala_scan_waterfall.png    # Hotspot residue waterfall chart
├── plot_ala_scan_heatmap.png      # ΔΔG heatmap
└── SUMMARY.txt                    # Full text report
```

---

## 8-Step Workflow

```
┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐
│ 1.Prepare├──►│2.Solvate ├──►│3.Minimize├──►│4.Equili- ├──►│5.Product-│
│  pdb2gmx │   │  +ions   │   │  steepest│   │  brate   │   │  ion MD  │
└──────────┘   └──────────┘   └──────────┘   │ NVT+NPT  │   │  (10 ns) │
                                              └──────────┘   └────┬─────┘
                                                                  │
                    ┌──────────┐   ┌──────────┐   ┌──────────────┘
                    │8.Analysis│◄──│7.Collect  │◄──│6. ASGBIE
                    │  plots   │   │  results  │   │  ├─ Phase 1: GB + IE
                    └──────────┘   └──────────┘   │  └─ Phase 2: Ala Scan
                                                  └──────────────────────
```

| Step | Script | Function | AMBER Equivalent |
|------|--------|----------|-----------------|
| 1 | `01_prepare.sh` | PDB → GROMACS topology | tLEaP |
| 2 | `02_solvate.sh` | Solvation + ions | tLEaP addions |
| 3 | `03_minimize.sh` | Energy minimization | pmemd min |
| 4 | `04_equilibrate.sh` | NVT + NPT equilibration | pmemd heat |
| 5 | `05_production.sh` | Production MD + RMSD | pmemd prod |
| 6 | `06_asgbie.sh` | **GB + IE + Alanine Scanning** | MMPBSA.py |
| 7 | `07_collect.sh` | Result collection | — |
| 8 | `08_analysis.py` | Plots & report | — |

### Run Options

```bash
bash run_asgbie.sh                 # Full workflow (Steps 1–8)
bash run_asgbie.sh --skip-md       # Skip MD, start from ASGBIE (Step 6)
bash run_asgbie.sh --only-asgbie   # Only ASGBIE calculation
bash run_asgbie.sh --step 4        # Resume from Step 4
```

---

## Batch Alanine Scanning

> **Key improvement over vanilla gmx_MMPBSA**: the built-in `&alanine_scanning` runs all residues in a single process — if one fails, everything fails. Our approach runs each residue independently.

### How It Works

```
For each residue in [ARG:25, GLU:30, TYR:45, ...]:
    1. Generate per-residue MMPBSA input file
    2. Run gmx_MMPBSA independently
    3. If already completed → skip (resume support)
    4. If failed → log error, continue to next
    5. Aggregate all ΔΔG → ala_scan_summary.csv
```

### Configuration

```bash
# config.sh

ALA_SCAN_ENABLED="true"                          # Toggle on/off
ALA_SCAN_RESIDUES="ARG:25,GLU:30,TYR:45,PHE:82"  # Manual list (recommended)
ALA_SCAN_RESIDUES=""                              # Empty = auto-detect all protein residues
ALA_SCAN_MAX_RESIDUES="50"                        # Cap for auto-detection
ALA_SCAN_PARALLEL="4"                             # Parallel jobs (requires GNU Parallel)
```

### Standalone Re-run

```bash
# Re-run failed residues (completed ones auto-skipped)
bash scripts/06b_ala_scan_batch.sh
```

### Output: `ala_scan_summary.csv`

| rank | residue | ddG_kcal_mol | hotspot |
|------|---------|-------------|---------|
| 1 | TYR:45 | +3.82 | YES |
| 2 | ARG:25 | +2.14 | YES |
| 3 | GLU:30 | +1.07 | YES |
| 4 | LEU:12 | +0.43 | no |

> Residues with **ΔΔG > 1.0 kcal/mol** are classified as hotspots.

---

## Configuration Reference

### Essential Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `INPUT_PDB` | complex.pdb | Input protein complex |
| `RECEPTOR_CHAINS` | A | Receptor chain ID(s) |
| `LIGAND_CHAINS` | B | Ligand chain ID(s) |
| `FORCEFIELD` | charmm27 | Force field (`charmm27` / `amber99sb-ildn` / `oplsaa`) |
| `PROD_STEPS` | 5000000 | Production MD steps (10 ns @ dt=0.002) |

### ASGBIE Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `IGB` | 5 | GB model (`5`=OBC2, recommended) |
| `IE_SEGMENT` | 25 | Interaction entropy segments |
| `REC_GROUP` | 1 | Receptor index group number ⚠️ |
| `LIG_GROUP` | 13 | Ligand index group number ⚠️ |

> ⚠️ **Important**: `REC_GROUP` and `LIG_GROUP` must match your system. Check with:
> ```bash
> gmx make_ndx -f asgbie_work/prod.tpr
> ```

### Simulation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `BOX_TYPE` | dodecahedron | Box shape (saves ~29% volume vs cubic) |
| `BOX_DIST` | 1.2 | Protein-to-edge distance (nm) |
| `ION_CONC` | 0.15 | NaCl concentration (M) |
| `TEMPERATURE` | 300 | Temperature (K) |
| `GMX_GPU` | "" | Set `"-nb gpu"` for GPU acceleration |
| `GMX_NTOMP` | 8 | OpenMP threads |

---

## Tips & Recommendations

### Choosing MD Length

| System Type | Recommended `PROD_STEPS` | Time |
|-------------|-------------------------|------|
| Small protein–peptide | 5,000,000 | 10 ns |
| Protein–protein | 25,000,000 | 50 ns |
| Antibody–antigen | 50,000,000 | 100 ns |

### Force Field Selection

| Force Field | Best For |
|-------------|----------|
| `charmm27` | General protein–protein (default) |
| `amber99sb-ildn` | Better backbone torsions, compatible with AMBER FF |
| `oplsaa` | Small molecule interactions |

### Performance Tips

- Use **GPU**: Set `GMX_GPU="-nb gpu"` in `config.sh`
- Use **parallel scanning**: Set `ALA_SCAN_PARALLEL` to half your CPU cores
- Use **dodecahedron box** (default): ~29% fewer water molecules than cubic
- **Reduce frames**: Set `MMPBSA_INTERVAL="2"` to analyze every 2nd frame

### Troubleshooting

| Issue | Solution |
|-------|---------|
| `gmx_MMPBSA` group number error | Run `gmx make_ndx -f prod.tpr`, update `REC_GROUP`/`LIG_GROUP` |
| Alanine scan residue failure | Check `ala_scan_logs/ala_*.log`, re-run `06b_ala_scan_batch.sh` |
| NVT temperature unstable | Increase `NVT_STEPS` to 100000 (200 ps) |
| High RMSD in production | System may need longer equilibration or position restraints |

---

## Project Structure

```
ASGBIE-GROMACS/
├── run_asgbie.sh              # One-click main script
├── config.sh                  # ← Edit this for your system
├── install_deps.sh            # Environment check / install
│
├── scripts/
│   ├── 01_prepare.sh          # PDB → GROMACS topology
│   ├── 02_solvate.sh          # Solvation + ions
│   ├── 03_minimize.sh         # Energy minimization
│   ├── 04_equilibrate.sh      # NVT + NPT equilibration
│   ├── 05_production.sh       # Production MD
│   ├── 06_asgbie.sh           # ASGBIE core (Phase 1: GB+IE, Phase 2: Ala scan)
│   ├── 06b_ala_scan_batch.sh  # Batch alanine scanning engine
│   ├── 07_collect.sh          # Result collection
│   └── 08_analysis.py         # Plotting & report generation
│
├── mdp/                       # GROMACS parameter files
│   ├── em.mdp                 # Energy minimization
│   ├── nvt.mdp                # NVT equilibration
│   ├── npt.mdp                # NPT equilibration
│   ├── prod.mdp               # Production MD
│   └── mmpbsa_asgbie.in       # MMPBSA input template
│
└── docs/
    └── REPORT.md              # Full technical report
```

## Software Dependencies

| Software | Version | Purpose | Install |
|----------|---------|---------|---------|
| [GROMACS](https://www.gromacs.org/) | ≥2022 | MD simulation | `conda install -c conda-forge gromacs` |
| [AmberTools](https://ambermd.org/AmberTools.php) | ≥22 | MMPBSA backend | `conda install -c conda-forge ambertools` |
| [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/) | ≥1.6 | GROMACS↔AMBER bridge | `pip install gmx_MMPBSA` |
| Python | ≥3.9 | Scripts & analysis | conda |
| [GNU Parallel](https://www.gnu.org/software/parallel/) | — | Parallel alanine scanning | `conda install -c conda-forge parallel` |

## References

1. Zhou M, et al. "An efficient approach to the accurate prediction of mutational effects in antigen binding to the MHC1." *Molecules* (2024) — **ASGBIE method**
2. Valdés-Tresanco MS, et al. "gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS." *JCTC* (2021) — **gmx_MMPBSA**
3. Li J, et al. "Significantly enhancing human antibody affinity via deep learning and ASGBIE." *Briefings in Bioinformatics* (2025) — **ASGBIE + AI**
4. Abraham MJ, et al. "GROMACS: High performance molecular simulations." *SoftwareX* (2015) — **GROMACS**
5. Original AMBER scripts: [xxttzhang/MD-ASGBIE](https://github.com/xxttzhang/MD-ASGBIE)

## License

[MIT](LICENSE)

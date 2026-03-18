#!/usr/bin/env python3
"""
ASGBIE-GROMACS Analysis & Plotting Script
Generates RMSD, energy decomposition, alanine scanning waterfall/heatmap plots
"""
import os
import sys
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ---------- Configuration ----------
RESULTS_DIR = os.environ.get("RESULTS_DIR", "asgbie_results")
WORKDIR = os.environ.get("WORKDIR", "asgbie_work")
DPI = 150
FIGSIZE = (10, 6)

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "figure.dpi": DPI,
})
sns.set_style("whitegrid")


def parse_xvg(filepath):
    """Parse a GROMACS .xvg file."""
    x, y = [], []
    if not os.path.exists(filepath):
        print(f"  [SKIP] {filepath} does not exist")
        return np.array([]), np.array([])
    with open(filepath) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
    return np.array(x), np.array(y)


def plot_rmsd():
    """Plot RMSD over time."""
    fig, ax = plt.subplots(figsize=FIGSIZE)

    for fname, label, color in [
        ("rmsd_backbone.xvg", "Backbone", "#2196F3"),
        ("rmsd_protein.xvg",  "Protein",  "#FF5722"),
    ]:
        t, rmsd = parse_xvg(os.path.join(RESULTS_DIR, fname))
        if len(t) > 0:
            ax.plot(t, rmsd, label=label, color=color, alpha=0.8, lw=1.5)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (nm)")
    ax.set_title("RMSD During Production MD")
    ax.legend()
    fig.tight_layout()
    outpath = os.path.join(RESULTS_DIR, "plot_rmsd.png")
    fig.savefig(outpath, dpi=DPI)
    plt.close()
    print(f"  [OK] {outpath}")


def plot_energy():
    """Plot energy minimization potential energy."""
    fig, ax = plt.subplots(figsize=FIGSIZE)
    t, pe = parse_xvg(os.path.join(RESULTS_DIR, "em_potential.xvg"))
    if len(t) > 0:
        ax.plot(t, pe, color="#4CAF50", lw=1.5)
        ax.set_xlabel("Step")
        ax.set_ylabel("Potential Energy (kJ/mol)")
        ax.set_title("Energy Minimization Convergence")
        fig.tight_layout()
        outpath = os.path.join(RESULTS_DIR, "plot_em_energy.png")
        fig.savefig(outpath, dpi=DPI)
        print(f"  [OK] {outpath}")
    plt.close()


def plot_equilibration():
    """Plot equilibration temperature and density."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    has_data = False

    # Temperature
    t, temp = parse_xvg(os.path.join(RESULTS_DIR, "nvt_temp.xvg"))
    if len(t) > 0:
        has_data = True
        axes[0].plot(t, temp, color="#E91E63", lw=0.8, alpha=0.7)
        axes[0].axhline(y=np.mean(temp), color="k", ls="--", lw=1)
        axes[0].set_xlabel("Time (ps)")
        axes[0].set_ylabel("Temperature (K)")
        axes[0].set_title(f"NVT Temperature (avg={np.mean(temp):.1f} K)")

    # Density
    t, dens = parse_xvg(os.path.join(RESULTS_DIR, "npt_density.xvg"))
    if len(t) > 0:
        has_data = True
        axes[1].plot(t, dens, color="#3F51B5", lw=0.8, alpha=0.7)
        axes[1].axhline(y=np.mean(dens), color="k", ls="--", lw=1)
        axes[1].set_xlabel("Time (ps)")
        axes[1].set_ylabel("Density (kg/m³)")
        axes[1].set_title(f"NPT Density (avg={np.mean(dens):.1f} kg/m³)")

    if has_data:
        fig.tight_layout()
        outpath = os.path.join(RESULTS_DIR, "plot_equilibration.png")
        fig.savefig(outpath, dpi=DPI)
        print(f"  [OK] {outpath}")
    plt.close()


def plot_asgbie_summary():
    """Parse and plot ASGBIE energy decomposition."""
    csv_path = os.path.join(RESULTS_DIR, "FINAL_RESULTS_ASGBIE.csv")
    if not os.path.exists(csv_path):
        csv_path = os.path.join(RESULTS_DIR, "FINAL_RESULTS_GBIE.csv")
    if not os.path.exists(csv_path):
        print(f"  [SKIP] ASGBIE CSV does not exist")
        return

    try:
        df = pd.read_csv(csv_path)
        print(f"  ASGBIE CSV: {len(df)} rows, columns: {list(df.columns[:8])}")

        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        if len(numeric_cols) > 0 and len(df) > 1:
            fig, ax = plt.subplots(figsize=(12, 6))
            energy_col = numeric_cols[0]
            ax.bar(range(len(df)), df[energy_col], color="#FF9800", alpha=0.8)
            ax.set_xlabel("Frame / Residue")
            ax.set_ylabel(f"{energy_col} (kcal/mol)")
            ax.set_title("ASGBIE Energy Decomposition")
            fig.tight_layout()
            outpath = os.path.join(RESULTS_DIR, "plot_asgbie_energy.png")
            fig.savefig(outpath, dpi=DPI)
            plt.close()
            print(f"  [OK] {outpath}")

    except Exception as e:
        print(f"  [WARN] Failed to parse CSV: {e}")


def plot_ala_scan_waterfall():
    """Plot alanine scanning waterfall chart (ddG ranking)."""
    # Find alanine scanning summary file
    ala_csv = os.path.join(WORKDIR, "asgbie_calc", "ala_scan_summary.csv")
    if not os.path.exists(ala_csv):
        ala_csv = os.path.join(RESULTS_DIR, "ala_scan_summary.csv")
    if not os.path.exists(ala_csv):
        print(f"  [SKIP] Alanine scanning results do not exist")
        return

    df = pd.read_csv(ala_csv)
    if len(df) == 0:
        print(f"  [SKIP] Alanine scanning results are empty")
        return

    print(f"  Alanine scanning: {len(df)} residues")

    # Sort by ddG descending, show Top 30
    df = df.sort_values("ddG_kcal_mol", ascending=False).head(30)

    # --- Waterfall plot ---
    fig, ax = plt.subplots(figsize=(14, 7))

    colors = []
    for ddg in df["ddG_kcal_mol"]:
        if ddg > 2.0:
            colors.append("#D32F2F")  # Strong hotspot (red)
        elif ddg > 1.0:
            colors.append("#FF5722")  # Medium hotspot (orange)
        elif ddg > 0:
            colors.append("#FFC107")  # Weak contribution (yellow)
        else:
            colors.append("#4CAF50")  # Unfavorable for binding (green)

    ax.bar(range(len(df)), df["ddG_kcal_mol"], color=colors, edgecolor="white", lw=0.5)

    # Add hotspot threshold lines
    ax.axhline(y=1.0, color="#FF5722", ls="--", lw=1, alpha=0.7, label="Hotspot threshold (1.0 kcal/mol)")
    ax.axhline(y=2.0, color="#D32F2F", ls="--", lw=1, alpha=0.7, label="Strong hotspot (2.0 kcal/mol)")
    ax.axhline(y=0, color="gray", ls="-", lw=0.5)

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["residue"], rotation=60, ha="right", fontsize=9)
    ax.set_xlabel("Residue")
    ax.set_ylabel("ΔΔG (kcal/mol)")
    ax.set_title("Alanine Scanning — Hotspot Residues (Top 30)")
    ax.legend(fontsize=10)

    fig.tight_layout()
    outpath = os.path.join(RESULTS_DIR, "plot_ala_scan_waterfall.png")
    fig.savefig(outpath, dpi=DPI)
    plt.close()
    print(f"  [OK] {outpath}")

    # --- Heatmap (if enough residues) ---
    if len(df) >= 3:
        width = min(30, max(10, len(df) * 0.5))
        fig, ax = plt.subplots(figsize=(width, 3))
        ddg_matrix = df[["ddG_kcal_mol"]].T
        ddg_matrix.columns = df["residue"].values

        sns.heatmap(
            ddg_matrix.astype(float),
            cmap="RdYlGn_r",
            center=0,
            annot=True,
            fmt=".1f",
            linewidths=0.5,
            ax=ax,
            cbar_kws={"label": "ΔΔG (kcal/mol)"},
        )
        ax.set_title("Alanine Scanning Heatmap")
        ax.set_ylabel("")
        fig.tight_layout()
        outpath = os.path.join(RESULTS_DIR, "plot_ala_scan_heatmap.png")
        fig.savefig(outpath, dpi=DPI)
        plt.close()
        print(f"  [OK] {outpath}")

    # Copy summary file to results directory
    dest = os.path.join(RESULTS_DIR, "ala_scan_summary.csv")
    if ala_csv != dest:
        shutil.copy(ala_csv, dest)


def generate_summary_report():
    """Generate text summary report."""
    report_path = os.path.join(RESULTS_DIR, "SUMMARY.txt")
    with open(report_path, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("  ASGBIE-GROMACS Analysis Report\n")
        f.write("=" * 60 + "\n\n")

        # RMSD statistics
        for fname in ["rmsd_backbone.xvg", "rmsd_protein.xvg"]:
            t, rmsd = parse_xvg(os.path.join(RESULTS_DIR, fname))
            if len(t) > 0:
                label = fname.replace("rmsd_", "").replace(".xvg", "")
                f.write(f"RMSD ({label}):\n")
                f.write(f"  Mean:   {np.mean(rmsd):.4f} nm\n")
                f.write(f"  StdDev: {np.std(rmsd):.4f} nm\n")
                f.write(f"  Max:    {np.max(rmsd):.4f} nm\n\n")

        # ASGBIE results
        for dat_name in ["FINAL_RESULTS_ASGBIE.dat", "FINAL_RESULTS_GBIE.dat"]:
            dat_path = os.path.join(RESULTS_DIR, dat_name)
            if os.path.exists(dat_path):
                f.write(f"GB + IE Results ({dat_name}):\n")
                with open(dat_path) as dat:
                    f.write(dat.read())
                f.write("\n")
                break

        # Alanine scanning results
        ala_csv = os.path.join(RESULTS_DIR, "ala_scan_summary.csv")
        if os.path.exists(ala_csv):
            df = pd.read_csv(ala_csv)
            hotspots = df[df["ddG_kcal_mol"].astype(float) > 1.0]
            f.write(f"\nAlanine Scanning:\n")
            f.write(f"  Total residues scanned: {len(df)}\n")
            f.write(f"  Hotspot residues (ddG > 1.0): {len(hotspots)}\n\n")
            if len(hotspots) > 0:
                f.write("  Hotspot Ranking:\n")
                f.write(f"  {'Rank':>4}  {'Residue':>12}  {'ddG':>10}\n")
                f.write(f"  {'-'*4}  {'-'*12}  {'-'*10}\n")
                for _, row in hotspots.iterrows():
                    f.write(f"  {int(row['rank']):>4}  {row['residue']:>12}  {float(row['ddG_kcal_mol']):>+10.2f}\n")
                f.write("\n")

        f.write("=" * 60 + "\n")
        f.write("Generated plots:\n")
        for png in sorted(os.listdir(RESULTS_DIR)):
            if png.endswith(".png"):
                f.write(f"  - {png}\n")

    print(f"  [OK] {report_path}")


# ========== Main ==========
if __name__ == "__main__":
    print("===== [Step 8] Analysis & Plotting =====")

    if not os.path.isdir(RESULTS_DIR):
        print(f"[ERROR] Results directory does not exist: {RESULTS_DIR}")
        sys.exit(1)

    plot_rmsd()
    plot_energy()
    plot_equilibration()
    plot_asgbie_summary()
    plot_ala_scan_waterfall()
    generate_summary_report()

    print("\n[Step 8] Done ✓")
    print(f"  All plots and reports saved to: {RESULTS_DIR}/")

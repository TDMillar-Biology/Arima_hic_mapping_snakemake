import os
import cooler
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Snakemake variables
mcool_path = snakemake.input.mcool
outdir = snakemake.params.outdir
sample_name = snakemake.params.sample_name

# Target chromosomes matching the directory structure
TARGET_CHROMS = ["2L", "2R", "3L", "3R", "4", "X", "Y"]

def plot_matrix(matrix, title, out_prefix):
    """Generates and saves a static matplotlib heatmap."""
    plt.figure(figsize=(10, 8))
    
    # Use LogNorm for visualization of absolute contact frequencies
    im = plt.imshow(
        matrix,
        cmap='viridis', # lets adjust this to pull from cli
        norm=LogNorm(vmax=np.nanmax(matrix)),
        interpolation='none'
    )
    
    plt.title(title)
    plt.colorbar(im, label='Absolute Contact Frequency')
    plt.xlabel("Bins")
    plt.ylabel("Bins")
    
    plt.tight_layout()
    plt.savefig(f"{out_prefix}.png", dpi=300)
    plt.savefig(f"{out_prefix}.pdf")
    plt.close()

def main():
    # 1. Plot Whole Genome
    # Use a coarser resolution for the whole genome to manage memory and image size
    # Assuming 100kb or 500kb exists in the mcool; fallback to base resolution if needed
    try:
        clr_wg = cooler.Cooler(f"{mcool_path}::/resolutions/100000")
    except KeyError:
        clr_wg = cooler.Cooler(f"{mcool_path}::/resolutions/10000")
    
    # Extract raw, absolute counts (balance=False)
    matrix_wg = clr_wg.matrix(balance=False, sparse=False)[:]
    wg_prefix = os.path.join(outdir, "hic_heatmap")
    plot_matrix(matrix_wg, f"{sample_name} - Whole Genome", wg_prefix)

    # 2. Plot Individual Chromosomes
    # Use the base resolution (10kb) for individual chromosomes
    clr_chrom = cooler.Cooler(f"{mcool_path}::/resolutions/10000")
    available_chroms = clr_chrom.chromnames

    for chrom in TARGET_CHROMS:
        if chrom in available_chroms:
            # Extract raw, absolute counts for the specific chromosome
            matrix_chrom = clr_chrom.matrix(balance=False, sparse=False).fetch(chrom)
            chrom_prefix = os.path.join(outdir, f"hic_heatmap_{chrom}")
            plot_matrix(matrix_chrom, f"{sample_name} - Chromosome {chrom}", chrom_prefix)

if __name__ == "__main__":
    main()
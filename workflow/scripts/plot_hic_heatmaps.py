import os
import copy
import cooler
import numpy as np
import matplotlib.pyplot as plt

# Snakemake variables mapped to Python
mcool_path = snakemake.input.mcool
outdir = snakemake.params.outdir
sample_name = snakemake.params.sample_name
cmap_name = snakemake.params.get("cmap", "viridis") 

# Target chromosomes matching the directory structure
TARGET_CHROMS = ["2L", "2R", "3L", "3R", "4", "X", "Y"]

# Dictionary defining the balancing states to extract and compare.
# 'False' yields the absolute raw counts. 'weight' applies the native ICE balancing.
BALANCING_METHODS = {
    "absolute": False, 
    "ICE_balanced": "weight"
}

def plot_matrix(matrix, title, out_prefix, extent, label):
    """Generates and saves a static matplotlib heatmap with genomic coordinates."""
    plt.figure(figsize=(10, 8))
    
    # 1. Apply the log1p transformation directly to the matrix array
    matrix_log = np.log1p(matrix)
    
    # 2. Calculate the 99th percentile for vmax on the TRANSFORMED matrix
    vmax = np.nanpercentile(matrix_log, 99)
    if np.isnan(vmax) or vmax <= 0:
        vmax = 1
        
    # 3. Safely load the colormap and color ICE-masked NaNs as the bottom colormap value
    try:
        current_cmap = copy.copy(plt.get_cmap(cmap_name))
    except ValueError:
        current_cmap = copy.copy(plt.get_cmap("viridis")) # Safe fallback
        
    # Set 'bad' (NaN) values to perfectly match the 0 value of the colormap
    current_cmap.set_bad(current_cmap(0.0))
        
    im = plt.imshow(
        matrix_log,
        cmap=current_cmap,
        vmax=vmax,
        vmin=0, # log1p(0) is 0
        interpolation='none',
        extent=extent # Maps matrix bins to genomic coordinates (Mbp)
    )
    
    plt.title(title)
    plt.colorbar(im, label=f"{label} (log1p)")
    plt.xlabel("Genomic Position (Mbp)")
    plt.ylabel("Genomic Position (Mbp)")
    
    plt.tight_layout()
    plt.savefig(f"{out_prefix}.png", dpi=300)
    plt.savefig(f"{out_prefix}.pdf")
    plt.close()

def main():
    # 1. Whole Genome Preparation
    try:
        clr_wg = cooler.Cooler(f"{mcool_path}::/resolutions/100000")
    except KeyError:
        clr_wg = cooler.Cooler(f"{mcool_path}::/resolutions/10000")
        
    # Calculate total genome length in Megabases for the plot extent
    wg_length_mb = clr_wg.chromsizes.sum() / 1e6
    wg_extent = [0, wg_length_mb, wg_length_mb, 0]

    # 2. Chromosome-Specific Preparation
    clr_chrom = cooler.Cooler(f"{mcool_path}::/resolutions/10000")
    available_chroms = clr_chrom.chromnames

    # 3. Iterate through balancing methods to generate comparisons
    for method_name, balance_col in BALANCING_METHODS.items():
        
        # Labeling for the colorbar
        cbar_label = 'Absolute Contact Frequency' if not balance_col else 'Balanced Contact Frequency'
        
        # --- Plot Whole Genome ---
        matrix_wg = clr_wg.matrix(balance=balance_col, sparse=False)[:]
        wg_prefix = os.path.join(outdir, f"hic_heatmap_WG_{method_name}")
        
        plot_matrix(
            matrix_wg, 
            f"{sample_name} - Whole Genome ({method_name})", 
            wg_prefix, 
            wg_extent, 
            cbar_label
        )

        # --- Plot Individual Chromosomes ---
        for chrom in TARGET_CHROMS:
            if chrom in available_chroms:
                matrix_chrom = clr_chrom.matrix(balance=balance_col, sparse=False).fetch(chrom)
                
                # Calculate single chromosome length in Megabases
                chrom_length_mb = clr_chrom.chromsizes[chrom] / 1e6
                chrom_extent = [0, chrom_length_mb, chrom_length_mb, 0]
                
                chrom_prefix = os.path.join(outdir, f"hic_heatmap_{chrom}_{method_name}")
                plot_matrix(
                    matrix_chrom, 
                    f"{sample_name} - Chromosome {chrom} ({method_name})", 
                    chrom_prefix, 
                    chrom_extent, 
                    cbar_label
                )

if __name__ == "__main__":
    main()
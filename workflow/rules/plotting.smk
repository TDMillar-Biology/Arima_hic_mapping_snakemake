# -----------------------------------------------------------------------------
# Step 10: Static Matplotlib Heatmap Generation
# -----------------------------------------------------------------------------
rule plot_hic_heatmaps:
    input:
        mcool = "results/{sample}/hic.mcool"
    output:
        wg_pdf = "results/{sample}/heatmaps/hic_heatmap_WG_absolute.pdf",
        wg_png = "results/{sample}/heatmaps/hic_heatmap_WG_absolute.png",
        wg_bal_pdf = "results/{sample}/heatmaps/hic_heatmap_WG_ICE_balanced.pdf"
    log:
        "logs/{sample}_plot_heatmaps.log"
    threads: 1
    conda: "../envs/cooler.yaml"
    params:
        outdir = "results/{sample}/heatmaps",
        sample_name = "{sample}",
        cmap = "RdYlBu_r" #viridis, magma, etc.
    script:
        "../scripts/plot_hic_heatmaps.py"
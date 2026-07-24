# -----------------------------------------------------------------------------
# Step 10: Static Matplotlib Heatmap Generation
# -----------------------------------------------------------------------------
rule plot_hic_heatmaps:
    input:
        mcool = "results/{sample}/hic.mcool"
    output:
        # Whole genome
        wg_pdf = "results/{sample}/hic_heatmap.pdf",
        wg_png = "results/{sample}/hic_heatmap.png",
        # Chromosome specific (using the expected Drosophila chromosomes as a dynamic flag)
        chrom_pdfs = expand("results/{{sample}}/hic_heatmap_{chrom}.pdf", chrom=["2L", "2R", "3L", "3R", "4", "X", "Y"]),
        chrom_pngs = expand("results/{{sample}}/hic_heatmap_{chrom}.png", chrom=["2L", "2R", "3L", "3R", "4", "X", "Y"])
    log:
        "logs/{sample}_plot_heatmaps.log"
    threads: 1
    conda: "../envs/cooler.yaml" # Assuming cooler env has matplotlib and numpy
    params:
        outdir = "results/{sample}",
        sample_name = "{sample}"
    script:
        "../scripts/plot_hic_heatmaps.py"
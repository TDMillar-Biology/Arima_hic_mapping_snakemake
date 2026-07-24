# -----------------------------------------------------------------------------
# Step 6: Extract Chromosome Sizes
# -----------------------------------------------------------------------------
rule get_chrom_sizes:
    input:
        ref_fai = config["reference"] + ".fai"
    output:
        sizes = "results/{sample}/chrom.sizes"
    shell:
        """
        cut -f1,2 {input.ref_fai} > {output.sizes}
        """

# -----------------------------------------------------------------------------
# Step 7: Parse BAM to Valid Pairs
# -----------------------------------------------------------------------------
# Since the BAM is coordinate-sorted and contains marked duplicates, we 
# stream it through samtools to drop duplicates (-F 1024), name-sort it on 
# the fly, and pipe the SAM output directly into pairtools.
rule bam_to_pairs:
    input:
        bam = "results/{sample}/hic_mapped.bam",
        chrom_sizes = "results/{sample}/chrom.sizes"
    output:
        pairs = "results/{sample}/hic.pairs"
    log:
        "logs/{sample}_bam_to_pairs.log"
    threads: config["threads"]
    conda: "../envs/pairtools.yaml"
    shell:
        """
        {config[samtools]} view -h -F 1024 -@ {threads} {input.bam} 2> {log} | \
        {config[samtools]} sort -n -@ {threads} -O SAM - 2>> {log} | \
        {config[pairtools]} parse -c {input.chrom_sizes} --drop-seq 2>> {log} | \
        {config[pairtools]} sort --nproc {threads} -o {output.pairs} 2>> {log}
        """

# -----------------------------------------------------------------------------
# Step 8: Build Base Cooler Matrix
# -----------------------------------------------------------------------------
rule cooler_cload:
    input:
        pairs = "results/{sample}/hic.pairs",
        chrom_sizes = "results/{sample}/chrom.sizes"
    output:
        cool = "results/{sample}/hic_10kb.cool"
    log:
        "logs/{sample}_cooler_cload.log"
    threads: config["threads"]
    conda: "../envs/cooler.yaml"
    shell:
        """
        {config[cooler]} cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
            {input.chrom_sizes}:{config[base_bin_size]} \
            {input.pairs} \
            {output.cool} &> {log}
        """

# -----------------------------------------------------------------------------
# Step 9: Generate Multi-Resolution Matrix
# -----------------------------------------------------------------------------
# The --no-balance flag is applied here to bypass iterative correction 
# algorithms (ICE/KR). This ensures the downstream matrices retain raw, 
# absolute contact frequencies rather than normalized weights.
rule cooler_zoomify:
    input:
        cool = "results/{sample}/hic_10kb.cool"
    output:
        mcool = "results/{sample}/hic.mcool"
    log:
        "logs/{sample}_cooler_zoomify.log"
    threads: config["threads"]
    conda: "../envs/cooler.yaml"
    shell:
        """
        {config[cooler]} zoomify \
            -n {threads} \
            -o {output.mcool} \
            {input.cool} &> {log}
        """
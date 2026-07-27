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
# Bypasses pairtools parse. Name-sorts the deduplicated Arima BAM, extracts 
# BEDPE format, and uses awk to print a standard 7-column pairs file for cooler.
rule bam_to_pairs:
    input:
        bam = "results/{sample}/hic_mapped.bam"
    output:
        pairs = "results/{sample}/hic.pairs"
    log:
        "logs/{sample}_bam_to_pairs.log"
    threads: config["threads"]
    resources:
        mem_mb = 16000
    conda: "../envs/matrix_prep.yaml"
    shell:
        """
        {config[samtools]} sort -n -@ {threads} {input.bam} 2> {log} | \
        bedtools bamtobed -bedpe -i - 2>> {log} | \
        awk -v OFS="\\t" '$1 != "." && $4 != "." {{print $7, $1, $2+1, $4, $5+1, $9, $10}}' > {output.pairs} 2>> {log}
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
        # Apply default ICE balancing natively
        {config[cooler]} zoomify \
            --balance \
            -n {threads} \
            -o {output.mcool} \
            {input.cool} &> {log}
        """
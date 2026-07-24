# -----------------------------------------------------------------------------
# Step 1: FASTQ to BAM
# -----------------------------------------------------------------------------
rule bwa_mem_r1:
    input:
        fq = lambda wildcards: samples_df.loc[wildcards.sample, "fq1"],
        ref = config["reference"]
    output:
        bam = temp("tmp/{sample}_1.bam")
    threads: config["threads"]
    conda: "../envs/hic_mapping.yaml"
    shell:
        """
        {config[bwa]} mem -t {threads} {input.ref} {input.fq} | \
        {config[samtools]} view -@ {threads} -Sb - > {output.bam}
        """

rule bwa_mem_r2:
    input:
        fq = lambda wildcards: samples_df.loc[wildcards.sample, "fq2"],
        ref = config["reference"]
    output:
        bam = temp("tmp/{sample}_2.bam")
    threads: config["threads"]
    conda: "../envs/hic_mapping.yaml"
    shell:
        """
        {config[bwa]} mem -t {threads} {input.ref} {input.fq} | \
        {config[samtools]} view -@ {threads} -Sb - > {output.bam}
        """

# -----------------------------------------------------------------------------
# Step 2: Filter 5' end
# -----------------------------------------------------------------------------
rule filter_5end:
    input:
        bam = "tmp/{sample}_{read}.bam",
        script = config["filter_script"]
    output:
        bam = temp("tmp/{sample}_{read}_filt.bam")
    conda: "../envs/hic_mapping.yaml"
    shell:
        """
        {config[samtools]} view -h {input.bam} | \
        perl {input.script} | \
        {config[samtools]} view -Sb - > {output.bam}
        """

# -----------------------------------------------------------------------------
# Step 3A: Pair reads & mapping quality filter
# -----------------------------------------------------------------------------
rule pair_and_sort:
    input:
        bam1 = "tmp/{sample}_1_filt.bam",
        bam2 = "tmp/{sample}_2_filt.bam",
        ref_fai = config["reference"] + ".fai",
        script = config["combiner_script"]
    output:
        bam = temp("tmp/{sample}_paired.bam")
    threads: config["threads"]
    conda: "../envs/hic_mapping.yaml"
    shell:
        """
        perl {input.script} {input.bam1} {input.bam2} {config[samtools]} {config[mapq_filter]} | \
        {config[samtools]} view -bS -t {input.ref_fai} - | \
        {config[samtools]} sort -@ {threads} -o {output.bam} -
        """

# -----------------------------------------------------------------------------
# Step 3B: Add read group
# -----------------------------------------------------------------------------
rule add_read_group:
    input:
        bam = "tmp/{sample}_paired.bam"
    output:
        bam = temp("tmp/{sample}_rg.bam")
    params:
        tmp_dir = "tmp/{sample}_rg_tmp",
        memory = config["rg_memory"]
    conda: "../envs/hic_mapping.yaml"
    shell:
        """
        mkdir -p {params.tmp_dir}
        {config[picard]} AddOrReplaceReadGroups \
            -Xmx{params.memory} \
            -Djava.io.tmpdir={params.tmp_dir} \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            ID={wildcards.sample} \
            LB={wildcards.sample} \
            SM={wildcards.sample} \
            PL=ILLUMINA \
            PU=none
        rm -rf {params.tmp_dir}
        """

# -----------------------------------------------------------------------------
# Step 4: Mark duplicates
# -----------------------------------------------------------------------------
rule mark_duplicates:
    input:
        bam = "tmp/{sample}_rg.bam"
    output:
        bam = "results/{sample}/hic_mapped.bam",
        metrics = "results/{sample}/metrics.{sample}.txt"
    params:
        tmp_dir = "tmp/{sample}_markdup_tmp",
        memory = config["markdup_memory"]
    conda: "../envs/hic_mapping.yaml"
    shell:
        """
        mkdir -p {params.tmp_dir}
        {config[picard]} MarkDuplicates \
            -Xmx{params.memory} \
            -XX:-UseGCOverheadLimit \
            -Djava.io.tmpdir={params.tmp_dir} \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            TMP_DIR={params.tmp_dir} \
            ASSUME_SORTED=TRUE \
            VALIDATION_STRINGENCY=LENIENT \
            REMOVE_DUPLICATES=TRUE
        rm -rf {params.tmp_dir}
        """

# -----------------------------------------------------------------------------
# Step 5: Index BAM and get stats
# -----------------------------------------------------------------------------
rule index_and_stats:
    input:
        bam = "results/{sample}/hic_mapped.bam",
        script = config["stats_script"]
    output:
        bai = "results/{sample}/hic_mapped.bam.bai",
        stats = "results/{sample}/hic_stats.txt"
    conda: "../envs/hic_mapping.yaml"
    shell:
        """
        {config[samtools]} index {input.bam}
        perl {input.script} {input.bam} > {output.stats}
        """
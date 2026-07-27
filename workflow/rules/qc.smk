# -----------------------------------------------------------------------------
# Step 11: Explicit Read Tracking Report
# -----------------------------------------------------------------------------
rule mapping_report:
    input:
        fq1 = lambda wildcards: samples_df.loc[wildcards.sample, "fq1"],
        bam = "results/{sample}/hic_mapped.bam"
    output:
        report = "results/{sample}/reports/mapping_report.txt"
    threads: 1
    conda: "../envs/matrix_prep.yaml" # Requires pysam
    script:
        "../scripts/mapping_report.py"

rule pipeline_attrition_report:
    input:
        fq1 = lambda wildcards: samples_df.loc[wildcards.sample, "fq1"],
        bwa1 = "tmp/{sample}_1.bam",
        bwa2 = "tmp/{sample}_2.bam",
        filt1 = "tmp/{sample}_1_filt.bam",
        filt2 = "tmp/{sample}_2_filt.bam",
        paired = "tmp/{sample}_paired.bam",
        final = "results/{sample}/hic_mapped.bam"
    output:
        report = "results/{sample}/reports/pipeline_attrition.tsv",
        plot_png = "results/{sample}/reports/pipeline_attrition.png",
        plot_pdf = "results/{sample}/reports/pipeline_attrition.pdf"
    threads: 4
    conda: "../envs/plotting.yaml" 
    script:
        "../scripts/plot_attrition.py"
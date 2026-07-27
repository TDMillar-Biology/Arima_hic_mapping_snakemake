import subprocess
import matplotlib.pyplot as plt

def get_count(command):
    """Executes a bash command and returns the integer output."""
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return int(result.stdout.strip())

def main():
    # Snakemake inputs
    fq1 = snakemake.input.fq1
    bwa1 = snakemake.input.bwa1
    bwa2 = snakemake.input.bwa2
    filt1 = snakemake.input.filt1
    filt2 = snakemake.input.filt2
    paired = snakemake.input.paired
    final = snakemake.input.final
    
    threads = snakemake.threads
    samtools = "samtools" # Assumes samtools is in the active conda env

    # 1. Absolute Raw Reads (R1 + R2)
    # wc -l divided by 4 gives read pairs; multiply by 2 for total reads
    raw_pairs = get_count(f"zcat {fq1} | wc -l") // 4
    raw_total = raw_pairs * 2

    # 2. BWA Mapped (Primary alignments only)
    bwa1_cnt = get_count(f"{samtools} view -c -@ {threads} -F 260 {bwa1}")
    bwa2_cnt = get_count(f"{samtools} view -c -@ {threads} -F 260 {bwa2}")
    bwa_total = bwa1_cnt + bwa2_cnt

    # 3. 5' Filtered
    filt1_cnt = get_count(f"{samtools} view -c -@ {threads} -F 260 {filt1}")
    filt2_cnt = get_count(f"{samtools} view -c -@ {threads} -F 260 {filt2}")
    filt_total = filt1_cnt + filt2_cnt

    # 4. Combiner Paired
    paired_total = get_count(f"{samtools} view -c -@ {threads} -F 260 {paired}")

    # 5. Final Deduplicated
    final_total = get_count(f"{samtools} view -c -@ {threads} -F 260 {final}")

    # Save TSV Report
    stages = ["Raw FASTQ", "BWA Mapped", "5' Filtered", "Combiner Paired", "Final Deduplicated"]
    counts = [raw_total, bwa_total, filt_total, paired_total, final_total]

    with open(snakemake.output.report, 'w') as f:
        f.write("Stage\tAbsolute_Reads\n")
        for stage, count in zip(stages, counts):
            f.write(f"{stage}\t{count}\n")

    # Generate Static Matplotlib Chart
    plt.figure(figsize=(10, 6))
    bars = plt.bar(stages, counts, color='black')
    
    plt.title(f"Hi-C Mapping Attrition - {snakemake.wildcards.sample}")
    plt.ylabel("Absolute Read Count")
    plt.xticks(rotation=45, ha="right")
    
    # Format y-axis to absolute integers (no scientific notation)
    plt.ticklabel_format(style='plain', axis='y')
    
    # Annotate absolute values atop bars
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval, f"{int(yval):,}", va='bottom', ha='center')

    plt.tight_layout()
    plt.savefig(snakemake.output.plot_png, dpi=300)
    plt.savefig(snakemake.output.plot_pdf)
    plt.close()

if __name__ == "__main__":
    main()
import gzip
import pysam
import sys

def main():
    fq1_path = snakemake.input.fq1
    bam_path = snakemake.input.bam
    out_path = snakemake.output.report
    sample = snakemake.wildcards.sample

    # 1. Establish absolute baseline from FASTQ
    # Every 4 lines is 1 read. Since this is paired-end, fq1 reads == total read pairs.
    raw_read_pairs = 0
    with gzip.open(fq1_path, 'rb') as f:
        for i, _ in enumerate(f):
            pass
        raw_read_pairs = (i + 1) // 4

    total_raw_reads = raw_read_pairs * 2

    # 2. Explicitly track hits in the BAM file
    # We are evaluating the BAM after duplicate removal.
    r1_mapped = 0
    r2_mapped = 0
    proper_pairs = 0
    interchromosomal_pairs = 0
    intrachromosomal_pairs = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            # Skip secondary or supplementary alignments to avoid double-counting
            if read.is_secondary or read.is_supplementary:
                continue

            if read.is_read1 and not read.is_unmapped:
                r1_mapped += 1
            if read.is_read2 and not read.is_unmapped:
                r2_mapped += 1

            if read.is_proper_pair:
                proper_pairs += 1
                
                # Check Hi-C specific mapping topology
                if read.reference_name == read.next_reference_name:
                    intrachromosomal_pairs += 1
                else:
                    interchromosomal_pairs += 1

    # Because Picard drops duplicates entirely from this BAM, these are absolute final counts.
    total_mapped_reads = r1_mapped + r2_mapped

    # 3. Write the explicit report
    with open(out_path, 'w') as out:
        out.write(f"Explicit Mapping Report for Sample: {sample}\n")
        out.write("="*50 + "\n")
        out.write(f"Absolute Raw Reads (FASTQ R1 + R2): {total_raw_reads}\n")
        out.write(f"Absolute Mapped Read 1 (BAM):     {r1_mapped}\n")
        out.write(f"Absolute Mapped Read 2 (BAM):     {r2_mapped}\n")
        out.write(f"Total Unique Mapped Reads:        {total_mapped_reads}\n")
        out.write("-" * 50 + "\n")
        out.write(f"Intra-chromosomal Contacts:       {intrachromosomal_pairs}\n")
        out.write(f"Inter-chromosomal Contacts:       {interchromosomal_pairs}\n")

if __name__ == "__main__":
    main()
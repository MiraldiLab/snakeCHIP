import os
from random import sample
import pandas as pd

configfile: "./inputs/config.yaml"

shell.prefix("set -eo pipefail; ")

# These can be put in the config.yaml file
OUT_DIR=config["output_dir"]
MEME_MOTIF_FILE=config["meme_motif_file"]
BLACKLIST=config["blacklist"]

class MetaTable(object):
    def __init__(self,
                 meta_path,
                 input_dir,
                 output_dir,
                 layout="SINGLE"):
        self.df = pd.read_table(meta_path)
        self.df = self.df[self.df["layout"] == layout]
        
        self.input_dir = input_dir
        self.output_dir = output_dir

        self.all_srr = self.df["srr"].unique().tolist()

        self.sample_list = self.df["srx"].unique().tolist()

    def getSRRList(self, wildcards):
        return self.df[self.df["srx"] == wildcards.sample]["srr"].tolist()

    def getReplicateFastq_SE(self, wildcards):
        SRR_LIST = self.df[self.df["srx"] == wildcards.sample]["srr"]
        
        return list(map(lambda srr_id: os.path.join(self.input_dir, "fasterq_dump", srr_id + ".fastq.gz"), SRR_LIST))

    def getReplicateFastq_pe1(self, wildcards):
        SRR_LIST = self.df[self.df["srx"] == wildcards.sample]["srr"]
        
        return list(map(lambda srr_id: os.path.join(self.input_dir, "fasterq_dump", srr_id + "_1.fastq.gz"), SRR_LIST))

    def getReplicateFastq_pe2(self, wildcards):
        SRR_LIST = self.df[self.df["srx"] == wildcards.sample]["srr"]
        
        return list(map(lambda srr_id: os.path.join(self.input_dir, "fasterq_dump", srr_id + "_2.fastq.gz"), SRR_LIST))

    def getTFName(self, wildcards):
        SRR_LIST = self.df[self.df["srx"] == wildcards.sample]["tf"]

meta = MetaTable(config["meta_file"], config["input_dir"], OUT_DIR)

# Expand the names out per sample

## Build the SRR file names for each replicate
SRR_SE = expand(config["input_dir"] + "/fasterq_dump/{accession}.fastq.gz", accession = meta.all_srr)

## Build the merged FASTQ file names
MERGED_FASTQ = expand(OUT_DIR + "/{sample}/merged_fastq/{sample}.fastq.gz", sample = meta.sample_list)

## Build the trimming FASTQ file names
TRIMMED_FASQ = expand(OUT_DIR + "/{sample}/trimmed_fastq/{sample}_trimmed.fq.gz", sample = meta.sample_list)

## Build names for SAM file
SAM_FILE = expand(OUT_DIR + "/{sample}/aligned_reads/{sample}.sam", sample = meta.sample_list)

## Build names for filtering and QC dependent on layout of library
NAMESORT_BAM_FILE = expand(OUT_DIR + "/{sample}/aligned_reads/{sample}.namesorted.bam", sample = meta.sample_list)

## Build names for deduplication and the final BAM files
ALL_DEDUPLICATED = expand(OUT_DIR + "/{sample}/aligned_reads/{sample}.deduplicated.bam", sample = meta.sample_list)
FINAL_BAM_FILE = expand(OUT_DIR + "/{sample}/aligned_reads/{sample}.bam", sample = meta.sample_list)
FINAL_POS_BAM_FILE = expand(os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_posStrand.bam"), sample = meta.sample_list)
FINAL_NEG_BAM_FILE = expand(os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_negStrand.bam"), sample = meta.sample_list)

# Optional Outputs
FILE_TYPES = ["postfiltered", "prefiltered", "final", "deduplicated","posStrandFinal", "negStrandFinal"]
ALL_FLAGSTAT = expand(OUT_DIR + "/{sample}/qc/flagstats/{sample}.bam.{types}.flagstat", sample = meta.sample_list, types = FILE_TYPES)
ALL_PEAKS = expand(OUT_DIR + "/{sample}/peaks/{sample}_ext147_peaks.narrowPeak", sample = meta.sample_list)
ALL_PEAKS_MACS_IDR = expand(OUT_DIR + "/{sample}/peaks/{sample}_ext147_p001_peaks_sorted.bed", sample = meta.sample_list)

ALL_TAGALIGN = expand(OUT_DIR + "/{sample}/tags/{sample}.bed.gz", sample = meta.sample_list)
FRIP = expand(os.path.join(OUT_DIR, "{sample}/qc/frip/{sample}_FRIP.txt"), sample = meta.sample_list)
FASTQC_posttrim=expand(os.path.join(OUT_DIR, "{sample}/qc/fastqc/post_trim/{sample}_trimmed.html"), sample = meta.sample_list)
ALL_HOMER = expand(os.path.join(OUT_DIR, "{sample}/homer/homerResults.html"), sample = meta.sample_list)
ALL_HOMER_SUMMARY = expand(os.path.join(OUT_DIR, "{sample}/homer/summary/{sample}_homer_summary.txt"), sample = meta.sample_list)
FINAL_BW_FILE = expand(os.path.join(OUT_DIR, "{sample}/bigwig/{sample}.bw"), sample = meta.sample_list)
FINAL_QUANT_RAW_BW_FILE = expand(os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_quantCPM.bw"), sample = meta.sample_list)
FINAL_POS_STRAND_QUANT_RAW_BW_FILE = expand(os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_posStrand_quantCPM.bw"), sample = meta.sample_list)
FINAL_NEG_STRAND_QUANT_RAW_BW_FILE = expand(os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_negStrand_quantCPM.bw"), sample = meta.sample_list)

FINAL_linearFE_BW_FILE = expand(os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_linearFE.bw"), sample = meta.sample_list)
FINAL_log10FE_BW_FILE = expand(os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_log10FE.bw"), sample = meta.sample_list)

FINAL_TSV = expand(os.path.join(OUT_DIR, "{sample}/hist/{sample}_chr22_arr_bp32_w32.tsv"), sample = meta.sample_list )
FINAL_SNS_TSV = expand(os.path.join(OUT_DIR, "{sample}/hist/{sample}_chr22_arr_bp32_w32_snsFormat.tsv"), sample = meta.sample_list)

FINAL_LINEARFE_TSV = expand(os.path.join(OUT_DIR, "{sample}/hist/{sample}_linearFE_chr22_arr_bp32_w32.tsv"), sample = meta.sample_list )
FINAL_LINEARFE_SNS_TSV = expand(os.path.join(OUT_DIR, "{sample}/hist/{sample}_linearFE_chr22_arr_bp32_w32_snsFormat.tsv"), sample = meta.sample_list )
FINAL_LOG10FE_TSV = expand(os.path.join(OUT_DIR, "{sample}/hist/{sample}_log10FE_chr22_arr_bp32_w32.tsv"), sample = meta.sample_list )
FINAL_LOG10FE_SNS_TSV = expand(os.path.join(OUT_DIR, "{sample}/hist/{sample}_log10FE_chr22_arr_bp32_w32_snsFormat.tsv"), sample = meta.sample_list )

rule all:
    input: FINAL_BW_FILE + ALL_TAGALIGN + FASTQC_posttrim + FRIP + ALL_FLAGSTAT + ALL_PEAKS + ALL_PEAKS_MACS_IDR + ALL_HOMER + FINAL_BAM_FILE + FINAL_POS_BAM_FILE + FINAL_NEG_BAM_FILE + FINAL_QUANT_RAW_BW_FILE + FINAL_POS_STRAND_QUANT_RAW_BW_FILE + FINAL_NEG_STRAND_QUANT_RAW_BW_FILE + FINAL_linearFE_BW_FILE + FINAL_log10FE_BW_FILE + FINAL_TSV + FINAL_SNS_TSV + FINAL_LINEARFE_TSV + FINAL_LINEARFE_SNS_TSV + FINAL_LOG10FE_TSV + FINAL_LOG10FE_SNS_TSV

rule get_fastq_se:
    priority: 2
    output:
        # the wildcard name must be accession, pointing to an SRA number
        os.path.join(config["input_dir"], "fasterq_dump/{accession}.fastq.gz")
    log:
        os.path.join(config["input_dir"], "logs/fasterq_dump/{accession}.log")
    threads: 8  # defaults to 6
    wrapper:
        "0.77.0/bio/sra-tools/fasterq-dump"

rule merge_replicates_fastq_SE:
    input: 
        meta.getReplicateFastq_SE

    output:
        temp(os.path.join(OUT_DIR, "{sample}/merged_fastq/{sample}.fastq.gz"))

    threads: 8
    message: "merging fastq files for R1: {input}"
    shell:
        """
        cat {input} > {output}
        """

rule fastqc_pre_trim:
    input:
        os.path.join(OUT_DIR, "{sample}/merged_fastq/{sample}.fastq.gz")
    output:
        html=os.path.join(OUT_DIR, "{sample}/qc/fastqc/pre_trim/{sample}.html"),
        zip=os.path.join(OUT_DIR, "{sample}/qc/fastqc/pre_trim/{sample}_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        os.path.join(OUT_DIR, "{sample}/logs/fastqc/{sample}_pre_trim.log")
    threads: 4
    wrapper:
        "0.77.0/bio/fastqc"

rule trim_galore_se:
    input:  
        os.path.join(OUT_DIR, "{sample}/merged_fastq/{sample}.fastq.gz")

    output:
        temp(os.path.join(OUT_DIR, "{sample}/trimmed_fastq/{sample}_trimmed.fq.gz"))

    threads: 6
    params: TRIM_DIR = os.path.join(OUT_DIR, "{sample}/trimmed_fastq")
    message: "trimming adaptors for {input} using trim_galore"
    log: os.path.join(OUT_DIR, "{sample}/logs/trimmed_fastq/{sample}.log")
    conda: "./envs/trim_galore.yaml"
    benchmark: os.path.join(OUT_DIR, "{sample}/logs/benchmark/{sample}_trim_adapter.benchmark")
    shell:
        """
        trim_galore -j {threads} -o {params.TRIM_DIR} {input} 2> {log}
        """

rule fastqc_post_trim:
    input:
        os.path.join(OUT_DIR, "{sample}/trimmed_fastq/{sample}_trimmed.fq.gz")
    output:
        html=os.path.join(OUT_DIR, "{sample}/qc/fastqc/post_trim/{sample}_trimmed.html"),
        zip=os.path.join(OUT_DIR, "{sample}/qc/fastqc/post_trim/{sample}_trimmed_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        os.path.join(OUT_DIR, "{sample}/logs/fastqc/{sample}_post_trim.log")
    threads: 4
    wrapper:
        "0.77.0/bio/fastqc"

rule bowtie2_align_SE:
    input:
        os.path.join(OUT_DIR, "{sample}/trimmed_fastq/{sample}_trimmed.fq.gz")

    output: temp(os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.sam"))
    
    threads: 8
    params: bowtie = "--very-sensitive --maxins 2000"
    message: "bowtie2_align_PE {input}: {threads} threads"
    conda: "./envs/bowtie2.yaml"
    benchmark: os.path.join(OUT_DIR, "{sample}/logs/benchmark/{sample}_bowtie2.benchmark")
    log: os.path.join(OUT_DIR, "{sample}/logs/bowtie2/{sample}.log")
    shell:
        """
        bowtie2 {params.bowtie} -p {threads} -x {config[idx_bt2]} -U {input} -S {output} 2> {log}
        """

# Get the flagstats for the alignment file
rule flagstat_prefiltered_bam:
    input:  os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.sam")

    output: os.path.join(OUT_DIR, "{sample}/qc/flagstats/{sample}.bam.prefiltered.flagstat")

    conda: "./envs/samtools.yaml"
    log:    os.path.join(OUT_DIR, "{sample}/logs/flagstat/{sample}.prefiltered.flagstat_bam")
    threads: 2
    message: "flagstat_prefiltered_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

rule quality_filter_namesort_sam2bam_se:
    input:  os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.sam")

    output: temp(os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.namesorted.bam"))

    log:    os.path.join(OUT_DIR, "{sample}/logs/filter_namesort/{sample}.namesort_bam")
    conda: "./envs/samtools.yaml"
    threads: 8
    benchmark: os.path.join(OUT_DIR, "{sample}/logs/benchmark/{sample}_filter_namesort.benchmark")
    message: "quality_filter_namesort_sam2bam {input}: {threads} threads"
    shell:
        """
        samtools view -@ {threads} -F 1804 -q 30 -b -u {input} | samtools sort -@ {threads} -n -o {output} - 2> {log}
        """

# Get the flagstats for the postfiltered BAM file
rule flagstat_postfiltered_bam:
    input:  os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.namesorted.bam")

    output: os.path.join(OUT_DIR, "{sample}/qc/flagstats/{sample}.bam.postfiltered.flagstat")

    log:    os.path.join(OUT_DIR, "{sample}/logs/flagstat/{sample}.postfiltered.flagstat_bam")
    conda: "./envs/samtools.yaml"
    threads: 2
    message: "flagstat_postfiltered_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

# Deduplicate the PCR duplicates in the BAM file using samtools markdup and samtools fixmate
rule deduplicate_BAM:
    input: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.namesorted.bam")

    output: temp(os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.fixmate.bam")), 
            temp(os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.deduplicated.bam")), 
            os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam"), 
            temp(os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.deduplicated.bam.bai")), 
            temp(os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.fixmate.bam.bai"))
    log:    os.path.join(OUT_DIR, "{sample}/logs/deduplicate/{sample}.fixmate_bam")
    threads: 8
    conda: "./envs/samtools.yaml"
    message: "deduplicate_BAM {input}: {threads} threads"
    shell:
        """
        samtools fixmate -@ {threads} -r -m {input} - 2> {log} | \
        samtools sort -@ {threads} -o {output[0]} - 2> {log}

        samtools index -@ {threads} {output[0]} 2> {log}

        samtools markdup -@ {threads} -r -s {output[0]} - 2> {log} | \
        samtools sort -@ {threads} -o {output[1]} - 2> {log}

        samtools index -@ {threads} {output[1]} 2> {log}

        samtools view -@ {threads} -b {output[1]} {config[keepChr]} | \
        samtools sort -@ {threads} -o {output[2]} - 2> {log}
        """

# Create Positive and Negative Strand
rule create_strand_BAM:
    input:
        os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam")
    output:
        os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_posStrand.bam"),
        os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_negStrand.bam")
    threads:
        4
    conda:
        "./envs/samtools.yaml"
    message:
        "Create Positive and Negative Strand BAM"
    shell:
        """
        samtools view -F 16 {input} -o {output[0]}
        samtools view -f 16 {input} -o {output[1]}
        """


# Create index 
rule index_final_BAM:
    input: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam")
    output: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam.bai")
    threads: 4
    conda: "./envs/samtools.yaml"
    message: "index final BAM {input}: {threads} threads"
    shell:
        """
        samtools index -@ {threads} {input}
        """

# Create index Positive Strand
rule index_final_positive_strand_BAM:
    input: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_posStrand.bam")
    output: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_posStrand.bam.bai")
    threads: 4
    conda: "./envs/samtools.yaml"
    message: "index positive strand final BAM {input}: {threads} threads"
    shell:
        """
        samtools index -@ {threads} {input}
        """

# Create index Negative Strand
rule index_final_negative_strand_BAM:
    input: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_negStrand.bam")
    output: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_negStrand.bam.bai")
    threads: 4
    conda: "./envs/samtools.yaml"
    message: "index negative strand final BAM {input}: {threads} threads"
    shell:
        """
        samtools index -@ {threads} {input}
        """

# Get the flagstats for the final bam
rule flagstat_deduplication_bam:
    input:  os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.deduplicated.bam")
    output: os.path.join(OUT_DIR, "{sample}/qc/flagstats/{sample}.bam.deduplicated.flagstat")
    log:    os.path.join(OUT_DIR, "{sample}/logs/flagstat/{sample}.deduplicated.flagstat_bam")
    threads: 2
    conda: "./envs/samtools.yaml"
    message: "flagstat_deduplicated_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

# Get the flagstats for the final bam
rule flagstat_final_bam:
    input:  os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam"),
            os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_posStrand.bam"),
            os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_negStrand.bam")

    output: os.path.join(OUT_DIR, "{sample}/qc/flagstats/{sample}.bam.final.flagstat"),
            os.path.join(OUT_DIR, "{sample}/qc/flagstats/{sample}.bam.posStrandFinal.flagstat"),
            os.path.join(OUT_DIR, "{sample}/qc/flagstats/{sample}.bam.negStrandFinal.flagstat")

    log:    os.path.join(OUT_DIR, "{sample}/logs/flagstat/{sample}.final.flagstat_bam"),
            os.path.join(OUT_DIR, "{sample}/logs/flagstat/{sample}.final.flagstat_posStrandFinal_bam"),
            os.path.join(OUT_DIR, "{sample}/logs/flagstat/{sample}.final.flagstat_negStrandFinal_bam")
    threads: 2
    conda: "./envs/samtools.yaml"
    message: "flagstat_final_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input[0]} > {output[0]} 2> {log[0]}

        samtools flagstat {input[1]} > {output[1]} 2> {log[1]}

        samtools flagstat {input[2]} > {output[2]} 2> {log[2]}
        """

# Convert alignments to bed file and bedpe file
rule convert_to_bed:
    input: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam")
    output: os.path.join(OUT_DIR, "{sample}/tags/{sample}.bed.gz")
    log:    os.path.join(OUT_DIR, "{sample}/logs/bamtobed/{sample}.tagalign")
    threads: 4
    conda: "./envs/bedtools.yaml"
    message: "tag_align {input}: {threads} threads"
    shell:
        """
        bedtools bamtobed -i {input} | gzip -nc > {output} 2> {log}
        """

# Call Peaks with MACS2
rule macs2_call_peaks:
    input: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam")

    output: os.path.join(OUT_DIR, "{sample}/peaks/{sample}_ext147_peaks.narrowPeak")

    log:    os.path.join(OUT_DIR, "{sample}/logs/macs2/{sample}.macs2")
    params: PEAK_DIR = os.path.join(OUT_DIR, "{sample}/peaks"),
            NAME = "{sample}_ext147"
    threads: 4
    conda: "./envs/macs2.yaml"
    message: "call peaks {input}: {threads} threads"
    shell:
        """
        macs2 callpeak -t {input} --name {params.NAME} -g hs --outdir {params.PEAK_DIR} --nomodel --extsize 147
        """

rule calculate_frip:
    input: 
           os.path.join(OUT_DIR, "{sample}/tags/{sample}.bed.gz"),
           os.path.join(OUT_DIR, "{sample}/peaks/{sample}_ext147_peaks.narrowPeak")

    output: os.path.join(OUT_DIR, "{sample}/qc/frip/{sample}_FRIP.txt")

    log:    os.path.join(OUT_DIR, "{sample}/logs/frip/{sample}.frip")
    threads: 4
    message: "calulate frip"
    conda: "./envs/bedtools.yaml"
    shell:
        """
        sh ./scripts/frip.sh {input[0]} {input[1]} {output}
        """

rule homer_analysis:
    input: PEAK_FILE = os.path.join(OUT_DIR, "{sample}/peaks/{sample}_ext147_peaks.narrowPeak"),
           MOTIF_FILE = MEME_MOTIF_FILE

    output: os.path.join(OUT_DIR, "{sample}/homer/homerResults.html")

    params: output_directory = os.path.join(OUT_DIR, "{sample}/homer")
    log:    os.path.join(OUT_DIR, "{sample}/logs/homer/{sample}.homer")
    threads: 4
    conda: "./envs/homer.yaml"
    message: "run homer motif analysis"
    shell:
        """
        # Run Homer for enrichment for known scanning
        findMotifsGenome.pl {input.PEAK_FILE} hg38 {params.output_directory} -mknown {input.MOTIF_FILE} -mcheck {input.MOTIF_FILE} -p 4

        # Run Homer for enrichemtn for de novo scanning
        #findMotifsGenome.pl {input.PEAK_FILE} hg38 {params.output_directory} -mcheck {input.MOTIF_FILE} -p 4
        """

# Call Peaks with MACS2
rule macs2_call_peaks_for_IDR:
    input: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam")

    output: temp(os.path.join(OUT_DIR, "{sample}/peaks/{sample}_ext147_p001_peaks.narrowPeak")),
            os.path.join(OUT_DIR, "{sample}/peaks/{sample}_ext147_p001_peaks_sorted.bed")

    log:    os.path.join(OUT_DIR, "{sample}/logs/macs2/{sample}.macs2_IDR")
    params: PEAK_DIR = os.path.join(OUT_DIR, "{sample}/peaks"),
            NAME = "{sample}_ext147_p001"
    threads: 4
    conda: "./envs/macs2.yaml"
    message: "call peaks {input}: {threads} threads"
    shell:
        """
        macs2 callpeak -t {input} --name {params.NAME} -g hs --outdir {params.PEAK_DIR} --nomodel --extsize 147 -B -p 1e-3

        sort -k8,8nr {output[0]} > {output[1]} 
        """

# Call Peaks with MACS2
rule macs2_call_peaks_for_IDR_posStrand:
    input: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_posStrand.bam")

    output: temp(os.path.join(OUT_DIR, "{sample}/peaks/{sample}_ext147_p001_peaks_posStrand.narrowPeak")),
            os.path.join(OUT_DIR, "{sample}/peaks/{sample}_ext147_p001_peaks_sorted_posStrand.bed")

    log:    os.path.join(OUT_DIR, "{sample}/logs/macs2/{sample}_posStrand.macs2_IDR")
    params: PEAK_DIR = os.path.join(OUT_DIR, "{sample}/peaks"),
            NAME = "{sample}_ext147_p001_posStrand"
    threads: 4
    conda: "./envs/macs2.yaml"
    message: "call peaks {input}: {threads} threads"
    shell:
        """
        macs2 callpeak -t {input} --name {params.NAME} -g hs --outdir {params.PEAK_DIR} --nomodel --extsize 147 -B -p 1e-3

        sort -k8,8nr {output[0]} > {output[1]}
        """


rule convert_bam_to_bigwig:
    input: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam"),
           os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam.bai")
    output: os.path.join(OUT_DIR, "{sample}/bigwig/{sample}.bw")
    log: os.path.join(OUT_DIR, "{sample}/logs/bigwig/{sample}.bigwig")
    threads: 4
    params: BLACKLIST
    conda: "./envs/deeptools.yaml"
    message: "convert {input} to bigwig: {threads} threads"
    shell:
        """
        bamCoverage --bam {input[0]} -o {output} --binSize 1 --normalizeUsing CPM -p {threads} --exactScaling -bl {params}
        """

rule convert_bam_to_bigwig_quantCPM:
    input: os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam"),
           os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}.bam.bai")
    output: os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_quantCPM.bw")
    log: os.path.join(OUT_DIR, "{sample}/logs/bigwig/{sample}.bigwig")
    threads: 4
    params: BLACKLIST
    conda: "./envs/deeptools.yaml"
    message: "convert {input} to bigwig: {threads} threads"
    shell:
        """
        bamCoverage --bam {input[0]} -o {output[0]} --binSize 1 --normalizeUsing CPM -p {threads} --exactScaling -bl {params}
        """

rule convert_bam_to_bigwig_posStrand_quantCPM:
    input:  os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_posStrand.bam"),
            os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_posStrand.bam.bai")
    output: os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_posStrand_quantCPM.bw")
    log: os.path.join(OUT_DIR, "{sample}/logs/bigwig/{sample}_posStrand.bigwig")
    threads: 4
    params: BLACKLIST
    conda: "./envs/deeptools.yaml"
    message: "convert {input} to bigwig: {threads} threads"
    shell:
        """
        bamCoverage --bam {input[0]} -o {output} --binSize 1 --normalizeUsing CPM -p {threads} --exactScaling -bl {params}
        """
        
rule convert_bam_to_bigwig_negStrand_quantCPM:
    input:  os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_negStrand.bam"),
            os.path.join(OUT_DIR, "{sample}/aligned_reads/{sample}_negStrand.bam.bai")
    output: os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_negStrand_quantCPM.bw")
    log: os.path.join(OUT_DIR, "{sample}/logs/bigwig/{sample}_negStrand.bigwig")
    threads: 4
    params: BLACKLIST
    conda: "./envs/deeptools.yaml"
    message: "convert {input} to bigwig: {threads} threads"
    shell:
        """
        bamCoverage --bam {input[0]} -o {output} --binSize 1 --normalizeUsing CPM -p {threads} --exactScaling -bl {params}
        """


rule make_bw_hist_file:
    input:  os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_quantCPM.bw")
    output: os.path.join(OUT_DIR, "{sample}/hist/{sample}_chr22_arr_bp32_w32.tsv"),
            os.path.join(OUT_DIR, "{sample}/hist/{sample}_chr22_arr_bp32_w32_snsFormat.tsv")
    log:    os.path.join(OUT_DIR, "{sample}/logs/hist/{sample}.makeTSV")
    threads: 4
    conda: "./envs/bwTOsnstsv.yaml"
    message: "Creating SNS format DF for histogram"
    shell:
        """
        python ./scripts/bwTOsnstsv.py --input_bw {input[0]} --out_name {output[0]} {output[1]}
        """
     

rule macs2_bdgcmp:
    input:  os.path.join(OUT_DIR, "{sample}/peaks/{sample}_ext147_p001_control_lambda.bdg"),
            os.path.join(OUT_DIR, "{sample}/peaks/{sample}_ext147_p001_treat_pileup.bdg")
    output: os.path.join(OUT_DIR, "{sample}/peaks/{sample}_linearFE.bdg"),
            os.path.join(OUT_DIR, "{sample}/peaks/{sample}_log10FE.bdg")
    log:    os.path.join(OUT_DIR, "{sample}/logs/macs2/{sample}.macs2_bdgcmp")
    threads: 4
    conda: "./envs/macs2_bdgcmp.yaml"
    message: "Create bedgraph with macs2 compare bedgraphs"
    shell:
        """
        macs2 bdgcmp -t {input[1]} -c {input[0]} -o {output[0]} -m FE -p 1
        macs2 bdgcmp -t {input[1]} -c {input[0]} -o {output[1]} -m logFE -p 1
        """


rule slop_bedGraph:
    input:  os.path.join(OUT_DIR, "{sample}/peaks/{sample}_linearFE.bdg"),
            os.path.join(OUT_DIR, "{sample}/peaks/{sample}_log10FE.bdg")
    output: os.path.join(OUT_DIR, "{sample}/peaks/{sample}_linearFE.fc.signal.bedgraph"),
            os.path.join(OUT_DIR, "{sample}/peaks/{sample}_log10FE.fc.signal.bedgraph")
    log: os.path.join(OUT_DIR, "{sample}/logs/peaks/{sample}_signal.bedgraph")
    threads: 4
    params: chrm_sizes = "/fs/ess/PES0738/20220614_maxatac_v1_data/snakemake/chip/inputs/hg38.chrom.sizes"
    conda: "./envs/slopBed.yaml"
    message: "macs bedgraph compare clean up step"
    shell:
        """
        slopBed -i {input[0]} -g {params.chrm_sizes} -b 0 | bedClip stdin {params.chrm_sizes} {output[0]}
        slopBed -i {input[1]} -g {params.chrm_sizes} -b 0 | bedClip stdin {params.chrm_sizes} {output[1]}

        rm -f {input[0]}
        rm -f {input[1]}
        """


rule sort_bedgraph:
    input:  os.path.join(OUT_DIR, "{sample}/peaks/{sample}_linearFE.fc.signal.bedgraph"),
            os.path.join(OUT_DIR, "{sample}/peaks/{sample}_log10FE.fc.signal.bedgraph")
    output: os.path.join(OUT_DIR, "{sample}/peaks/{sample}_linearFE.fc.signal.srt.bedgraph"),
            os.path.join(OUT_DIR, "{sample}/peaks/{sample}_log10FE.fc.signal.srt.bedgraph")
    threads: 4
    message: "Sorting Bedgraphs"
    shell:
        """
        sort -k1,1 -k2,2n {input[0]} > {output[0]}
        sort -k1,1 -k2,2n {input[1]} > {output[1]}
        """

rule bedGraph_To_Bigwig:
    input:
        os.path.join(OUT_DIR, "{sample}/peaks/{sample}_linearFE.fc.signal.srt.bedgraph"),
        os.path.join(OUT_DIR, "{sample}/peaks/{sample}_log10FE.fc.signal.srt.bedgraph")
    output:
        os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_linearFE.bw"),
        os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_log10FE.bw")
    log:
        os.path.join(OUT_DIR, "{sample}/logs/bigwig/{sample}.bedGraphToBigWig")
    threads: 
        4
    message:
        "Generate bedGraphToBigWig"
    params:
        chrm_sizes = "/fs/ess/PES0738/20220614_maxatac_v1_data/snakemake/chip/inputs/hg38.chrom.sizes"
    conda:
        "./envs/bedGraphToBigWig.yaml"
    shell:
        """
        bedGraphToBigWig {input[0]} {params.chrm_sizes} {output[0]}
        bedGraphToBigWig {input[1]} {params.chrm_sizes} {output[1]}

        rm -f {input[0]}
        rm -f {input[1]}
        """

rule make_FE_bw_hist_file:
    input:  os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_linearFE.bw"),
            os.path.join(OUT_DIR, "{sample}/bigwig/{sample}_log10FE.bw")
    output: os.path.join(OUT_DIR, "{sample}/hist/{sample}_linearFE_chr22_arr_bp32_w32.tsv"),
            os.path.join(OUT_DIR, "{sample}/hist/{sample}_linearFE_chr22_arr_bp32_w32_snsFormat.tsv"),
            os.path.join(OUT_DIR, "{sample}/hist/{sample}_log10FE_chr22_arr_bp32_w32.tsv"),
            os.path.join(OUT_DIR, "{sample}/hist/{sample}_log10FE_chr22_arr_bp32_w32_snsFormat.tsv")

    log:    os.path.join(OUT_DIR, "{sample}/logs/hist/{sample}_linear_log.makeTSV")
    threads: 4
    conda: "./envs/bwTOsnstsv.yaml"
    message: "Creating SNS format DF for histogram"
    shell:
        """
        python ./scripts/bwTOsnstsv.py --input_bw {input[0]} --out_name {output[0]} {output[1]}
        python ./scripts/bwTOsnstsv.py --input_bw {input[1]} --out_name {output[2]} {output[3]}
        """
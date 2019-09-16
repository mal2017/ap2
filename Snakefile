from os.path import realpath
from os.path import split as pathsplit
import subprocess

# info
__author__ = "Matt Lawlor"


# set shell
shell.executable("/bin/bash")

#singularity: "docker://continuumio/miniconda3:4.4.10"


## default
# prealignment to bait seqs
# remove chrom m and save
# filter discordant reads
# filter unmapped reads
# convert to cram
# save my reference genomes somewhere safe!!
# rule to download from sra, or ena, or ftp, or box, or dropbox, or google drive
# rule to combine multiple fqs per end
# align with bt2 or bwa aln
# idr
# normalize peak scores
# choose best peaks from cohort
# quantify and store as hdf5 or something compressed
# cuts per covered base
# tss plots
# frips
# methods: https://www.biostars.org/p/220268/
## special subworkflows
# samtools stats
# multiqc
# regions of interest plots
# nucleosome calling
# footprint calling
# copy number
# se calling?
# bigwig generation
#

rule target:
    input:
        #expand("alignments/{s}.raw.sam",s=config.get("samples",None))
        #"index.sr.mmi"
        #expand("index_bt2.{x}",x=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]),
        expand("alignments_bt2/{s}.raw.{f}",s=config.get("samples",None),f=["sam","bam","cram"])

rule index_genome_bt2:
    input:
        config.get("genome",None)
    output:
        expand("index_bt2.{x}",x=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"])
    singularity:
        "docker://quay.io/biocontainers/bowtie2:2.3.5--py37he860b03_0"
    conda:
        "environments/bowtie2.yaml"
    shell:
        "bowtie2-build {input} index_bt2"

rule align_bt2:
    input:
        r1 = lambda wc: config["samples"][wc.s]["fastq"]["r1"],
        r2 = lambda wc: config["samples"][wc.s]["fastq"]["r2"],
        idx = rules.index_genome_bt2.output,
    output:
        "alignments_bt2/{s}.raw.sam"
    conda:
        "environments/bowtie2.yaml"
    threads:
        4
    shell:
        "bowtie2 --trim-to 3:30 -p {threads} --phred33 "
        "--very-fast-local -X 2000 "
        "-x index_bt2 -1 {input.r1} -2 {input.r2} "
        "-S {output}"


rule index_genome_mm2:
    input:
        config.get("genome",None)
    output:
        "index_mm2.sr.mmi"
    conda:
        "environments/minimap2.yaml"
    shell:
        "minimap2 -x sr -d {output} {input}"

rule align_mm2:
    input:
        r1 = lambda wc: config["samples"][wc.s]["fastq"]["r1"],
        r2 = lambda wc: config["samples"][wc.s]["fastq"]["r2"],
        idx = rules.index_genome_mm2.output,
    output:
        "alignments_mm2/{s}.raw.sam"
    conda:
        "environments/minimap2.yaml"
    threads:
        1
    shell:
        """minimap2 -ax sr {input.idx} -o {output} -t {threads} """
        """{input.r1} {input.r2}"""

rule sam_to_bam:
    input:
        sam="{file}.sam",
    output:
        "{file}.bam"
    threads:
        2
    shell:
        "samtools sort {input.sam} -@ {threads} -o {output}"

rule unzip_genome:
    input:
        g=config.get("genome",None)
    output:
        temp("temp_genome.fa")
    shell:
        "gzip --keep -d {input.g} -c > {output}"

rule sam_to_cram:
    input:
        sam="{file}.sam",
        g=rules.unzip_genome.output
    output:
        "{file}.cram"
    threads:
        2
    shell:
        "samtools sort {input.sam} -@ {threads} -o {output} --reference {input.g}"

rule cram_to_bam:
    input:
        cram="{file}.cram",
        g=rules.unzip_genome.output
    output:
        "{file}.bam"
    shell:
        "samtools view -b {input.cram} -o {output} -T {input.g}"

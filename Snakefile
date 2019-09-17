from os.path import realpath
from os.path import split as pathsplit
import subprocess

# info
__author__ = "Matt Lawlor"

# set shell
shell.executable("/bin/bash")

configfile: "config.yaml"

"""
CLUSTER_NAME=snk-cl2
NODES=2
ZONE=us-central1-a
REMOTE=GS
PREFIX=archibald
MACHINE_TYPE=n1-standard-4
gcloud container clusters create $CLUSTER_NAME \
    --num-nodes=$NODES \
    --scopes storage-rw \
    --machine-type=$MACHINE_TYPE \
    --zone $ZONE
gcloud container clusters get-credentials $CLUSTER_NAME --zone $ZONE

snakemake --kubernetes --use-conda \
    --default-remote-provider $REMOTE \
    --default-remote-prefix $PREFIX \
    --latency-wait 300

# after
gcloud container clusters delete $CLUSTER_NAME --zone $ZONE
"""


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
        expand("{s}.subs.cram",s=config.get("samples",None)),

#http://biolearnr.blogspot.com/2017/11/snakemake-using-inputoutput-values-in.html
# think about using this to substitute prefix so that the prefix includes archibald
# or whatever the remote should be.
rule index_genome_bt2:
    input:
        config.get("genome",None)
    output:
        expand("{g}.{x}",g= config.get("genome",None), x=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"])
    singularity:
        "docker://quay.io/biocontainers/bowtie2:2.3.5--py37he860b03_0"
    conda:
        "environments/bowtie2.yaml"
    shell:
        "bowtie2-build {input} {input}"

rule align_bt2:
    input:
        r1 = lambda wc: config["samples"][wc.s]["fastq"]["r1"],
        r2 = lambda wc: config["samples"][wc.s]["fastq"]["r2"],
        idx = rules.index_genome_bt2.output,
        idx_path = config.get("genome",None)
    output:
        "{s}.raw.sam"
    conda:
        "environments/bowtie2.yaml"
    threads:
        2
    shell:
        "bowtie2 --trim-to 3:30 -p {threads} --phred33 "
        "--very-fast-local -X 2000 "
        "-x {input.idx_path} -1 {input.r1} -2 {input.r2} "
        "-S {output}"

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
    conda:
        "environments/samtools.yaml"
    threads:
        2
    shell:
        "samtools sort {input.sam} -@ {threads} -o {output} --reference {input.g}"

rule cram_to_bam:
    input:
        cram="{file}.cram",
        g=rules.unzip_genome.output
    output:
        temp("{file}.bam")
    conda:
        "environments/samtools.yaml"
    shell:
        "samtools view -b {input.cram} -o {output} -T {input.g}"

rule subsamp:
    input:
        "{s}.raw.cram"
    output:
        "{s}.subs.cram"
    conda:
        "environments/samtools.yaml"
    shell:
        "samtools view -C {input} -o {output} -s 0.0005"

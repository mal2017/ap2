from os.path import realpath
from os.path import split as pathsplit
import subprocess
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

# info
__author__ = "Matt Lawlor"

# set shell
shell.executable("/bin/bash")

# configuration for cloud runs
configfile: "config.yaml"

# set remotes
GS = GSRemoteProvider()
FTP = FTPRemoteProvider()

# remote resources
hg38_idx_files = ["seq-resources/NCBI-hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index." + x for x in ["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]]
hg38_idx_paths = [GS.remote(y) for y in hg38_idx_files]
hg38_fa = GS.remote("seq-resources/NCBI-hg38/hg38.fa.gz")
hg38_fai = GS.remote("seq-resources/NCBI-hg38/hg38.fa.gz.fai")
hg38_gzi = GS.remote("seq-resources/NCBI-hg38/hg38.fa.gz.gzi")
hg38_idx_pfx = "seq-resources/NCBI-hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"

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


rule target:
    input:
        expand("{s}.raw.cram",s=config.get("samples",None)),


#http://biolearnr.blogspot.com/2017/11/snakemake-using-inputoutput-values-in.html
rule align_bt2:
    input:
        r1 = lambda wc: FTP.remote(config["samples"][wc.s]["fastq"]["r1"]),
        r2 = lambda wc: FTP.remote(config["samples"][wc.s]["fastq"]["r2"]),
        idx = hg38_idx_paths,
        fa=hg38_fa,
        fai=hg38_fai,
        gzi=hg38_gzi,
    output:
        temp("{s}.raw.cram")
    conda:
        "environments/bowtie2.yaml"
    threads:
        2
    singularity:
        "docker://quay.io/biocontainers/bowtie2:2.3.5--py37he860b03_0"
    params:
        idx_pfx = hg38_idx_pfx,
    shell:
        "bowtie2 --trim-to 3:30 --phred33 "
        "--no-discordant "
        "--no-unal "
        "-k 1 -X 800 "
        "-x {params.idx_pfx} -1 {input.r1} -2 {input.r2} "
        #"-S {output} "
        "--skip 1000000 -u 10000000 | "
        "samtools sort -O cram -o {output} --reference {input.fa}"

# rule sam_to_cram:
#     input:
#         sam="{file}.sam",
#         fa=hg38_fa,
#         fai=hg38_fai,
#         gzi=hg38_gzi,
#     output:
#         "{file}.cram"
#     conda:
#         "environments/samtools.yaml"
#     threads:
#         2
#     shell:
#         "samtools sort {input.sam} -@ {threads} -O cram -o {output} --reference {input.fa}"

rule cram_to_bam:
    input:
        cram="{file}.cram",
        fa=hg38_fa,
        fai=hg38_fai,
        gzi=hg38_gzi,
    output:
        temp("{file}.bam")
    conda:
        "environments/samtools.yaml"
    shell:
        "samtools view -b {input.cram} -o {output} -T {input.fa}"

rule subsamp:
    input:
        c="{s}.raw.cram",
        fa=hg38_fa,
        fai=hg38_fai,
        gzi=hg38_gzi,
    output:
        "{s}.subs.cram"
    conda:
        "environments/samtools.yaml"
    shell:
        "samtools view -C {input.c} -o {output} -s 0.0005 -T {input.fa}"

## TODO
# function for making remote objects or checking if local and either making a remote object or a path including sra
# prealignment to bait seqs
# filter unwanted chroms
# rule to combine multiple fqs per end
# rule to count kmers? motifs? maybe convert to fasta first?




# idr
# normalize peak scores
# choose best peaks from cohort
# quantify and store as sparse hdf5 or something compressed
# cuts per covered base per total reads
# histogram/distribution of cuts/base
# frips
# make a dataframe with metadata
# methods: https://www.biostars.org/p/220268/

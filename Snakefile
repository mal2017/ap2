from os.path import realpath
from os.path import split as pathsplit
import subprocess
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import sys

# Block annoying warnings
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# META
__author__ = "Matt Lawlor"

# SETUP
shell.executable("/bin/bash")
GS = GSRemoteProvider()
FTP = FTPRemoteProvider()
S3 = S3RemoteProvider()
HTTP = HTTPRemoteProvider()
configfile: "config.yaml"

# PARAMS
BT2_TRIM_SIZE = config.get("BT2_TRIM_SIZE",30)
BT2_MAX_ISIZE = config.get("BT2_MAX_ISIZE",800)
MACS2_SHIFT_SIZE = config.get("MACS2_SHIFT_SIZE",-100)
MACS2_EXTENSION = config.get("MACS2_EXTENSION",200)
PYDNASE_PVAL = config.get("PYDNASE_PVAL",-2)
MAPQ_CUTOFF = config.get("MAPQ_CUTOFF",30)
GENOME_SIZE = config.get("GENOME_SIZE",None)

# DETERMINE REMOTE OR LOCAL RESOURCE
def determine_resource(path):
    if "gs://" in path:
         return GS.remote(path.replace("gs://",""))
    elif "ftp://" in path:
         return FTP.remote(path)
    elif "s3://" in path:
         return S3.remote(path.replace("s3://",""))
    elif "http://" in path:
         return HTTP.remote(path.replace("http://",""))
    elif "https://" in path:
         return HTTP.remote(path.replace("https://",""))
    else:
        return path

# REMOTE RESOURCES
bt2_idx_paths = [determine_resource(y) for y in config.get("BT2_FILES",None)]
genome_fa = determine_resource(config.get("GENOME_FA",None))
genome_fai = determine_resource(config.get("GENOME_FAI",None))
genome_gzi = determine_resource(config.get("GENOME_GZI",None))
genome_chroi_names = config.get("GENOME_CHR",None).split(",")
genome_bl = determine_resource(config.get("GENOME_BL",None))


rule target:
    """
    Generate footprints, peaks, and a cut matrix from all samples in the manifest.
    """
    input:
        expand("{s}_peaks.mrg.bed",s=config.get("samples",None)),
        expand("{s}_smts.bed",s=config.get("samples",None)),
        expand("{s}_fps.bed",s=config.get("samples",None)),
        expand("{s}.cpm.bw",s=config.get("samples",None)),

# ------------------------------------------------------------------------------
# Preproc
# ------------------------------------------------------------------------------

#http://biolearnr.blogspot.com/2017/11/snakemake-using-inputoutput-values-in.html
rule align_bt2:
    """
    Align reads with bowtie2.
    """
    input:
        r1 = lambda wc: [determine_resource(x) for x in config["samples"][wc.s]["fastq"]["r1"]],
        r2 = lambda wc: [determine_resource(x) for x in config["samples"][wc.s]["fastq"]["r2"]],
        idx = bt2_idx_paths,
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
    output:
        temp("{s}.raw.cram")
    conda:
        "envs/bowtie2.yaml"
    threads:
        1
    params:
        idx_pfx = config.get("BT2_IDX_PFX",None),
        trim=BT2_TRIM_SIZE,
        isize=BT2_MAX_ISIZE
    shell:
        "bowtie2 --trim-to 3:{params.trim} --phred33 "
        "--no-discordant "
        "--no-unal "
        "-k 1 -X {params.isize} "
        "-x {params.idx_pfx} -1 {input.r1} -2 {input.r2} |"
        #"--skip 1000000 "
        #"-u 10000000 | "
        "samtools sort -O cram --output-fmt-option lossy_names=1,level=9,store_md=0,store_nm=0 "
        "-o {output} --reference {input.fa}"

rule blacklist_filter_reads:
    """
    Remove blacklistlisted reads.
    """
    input:
        crm=rules.align_bt2.output,
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
        bl=genome_bl,
        crai="{s}.raw.cram.crai"
    output:
        temp("{s}.bl.cram")
    conda:
        "envs/bedtools.yaml"
    threads:
        1
    shell:
        "CRAM_REFERENCE={input.fa} "
        "bedtools intersect -v -a {input.crm} -b {input.bl} > {output}"

rule fix_mate_info:
    """
    Update mate info in aux tags and coord sort.
    """
    input:
        crm="{s}.bl.nsrt.cram",
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
    output:
        temp("{s}.fixm.cram")
    conda:
        "envs/bowtie2.yaml"
    threads: 1
    shell:
        "samtools fixmate -m --reference {input.fa} {input.crm} - | "
        "samtools sort -O cram "
        "--output-fmt-option lossy_names=1,level=9,store_md=0,store_nm=0 "
        "-o {output} --reference {input.fa}"

rule clean_reads:
    """
    Filter reads by MAPQ, canonical chromosomes, pcr dups.
    """
    input:
        crm="{s}.fixm.cram",
        crai="{s}.fixm.cram.crai",
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
    output:
        "{s}.clean.cram"
    params:
        chr=genome_chroi_names,
        mapq=MAPQ_CUTOFF,
    conda:
        "envs/bowtie2.yaml"
    threads:
        1
    shell:
        "samtools view -u -q {params.mapq} "
        "--reference {input.fa} {input.crm} {params.chr} | "
        "samtools markdup -r "
        "--output-fmt-option lossy_names=1,level=9,store_md=0,store_nm=0 "
        "--reference {input.fa} - {output}"

# ------------------------------------------------------------------------------
# Peak calling
# ------------------------------------------------------------------------------

rule call_peaks:
    """
    Call peaks using macs2 and shifting the reads to center on the cut site.
    """
    input:
        crm="{s}.clean.bam",
    output:
        sum="{s}_smts.bed",
        #np="{s}_peaks.narrowPeak",
        mrg="{s}_peaks.mrg.bed"
    params:
        gs=GENOME_SIZE,
        shift=MACS2_SHIFT_SIZE,
        ext=MACS2_EXTENSION
    conda:
        "envs/macs2.yaml"
    shadow: "shallow"
    threads:
        1
    shell:
        "macs2 callpeak -t {input} -f BAM " # TODO right now this only takes the 1st read in the pair, which is find for now...
        "--nomodel --shift {params.shift} --keep-dup all "
        "--extsize {params.ext} "
        "-g {params.gs} -n {wildcards.s} --call-summits; "
        "bedtools sort -i {wildcards.s}_peaks.narrowPeak | bedtools merge > {output.mrg}; "
        "cp {wildcards.s}_summits.bed {output.sum}"

# ------------------------------------------------------------------------------
# Footprinting
# ------------------------------------------------------------------------------

rule call_footprints:
    """
    Call footprints.
    """
    input:
        bam="{s}.clean.bam",
        bai="{s}.clean.bam.bai",
        bed="{s}_peaks.mrg.bed"
    output:
        "{s}_fps.bed"
    shadow:
        "shallow"
    conda:
        "envs/pydnase.yaml"
    threads:
        1
    params:
        odir= "{s}-fp-tmp/",
        pv = PYDNASE_PVAL
    shell:
        "rm -rf {params.odir} && mkdir {params.odir} && "
        "wellington_footprints.py -p {threads} -o {wildcards.s} -pv {params.pv} -A {input.bed} {input.bam} {params.odir}; "
        "mv {params.odir}/p\ value\ cutoffs/{wildcards.s}.WellingtonFootprints.{params.pv}.bed {output}; "
        "rm -rf {params.odir}"

# ------------------------------------------------------------------------------
# Footprinting
# ------------------------------------------------------------------------------

rule make_bigwigs:
    input:
        crm="{s}.clean.cram",
        crai="{s}.clean.cram.crai"
    output:
        "{s}.cpm.bw"
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input.crm} "
        "--Offset 4 6 --outFileName {output} "
        "--outFileFormat bigwig "
        "--binSize 50 --smoothLength 150 "
        "--verbose --normalizeUsing CPM"

## TODO
# for sra, have a rule that pipes directly from sra-dump to bowtie2
# function for making remote objects or checking if local and either making a remote object or a path including sra
#       Might have to prioritize an sra dump rule first?
# rule to combine multiple fqs per end
# bamCoverage -b GM12878.clean.cram --Offset 4 6 --outFileName GM12878.cpm.bw --outFileFormat bigwig --binSize 50 --smoothLength 150 --verbose --normalizeUsing CPM
# idr or normalize peak scores like the corces atac paper

## DOWNSTREAM STUFF
# choose best peaks from cohort
# quantify and store as sparse hdf5 or something compressed
# cuts per covered base per total reads
# histogram/distribution of cuts/base
# frips
# make a dataframe with metadata
# methods: https://www.biostars.org/p/220268/

# ------------------------------------------------------------------------------
# Generics
# ------------------------------------------------------------------------------

rule index_cram:
    input:
        "{file}.cram"
    output:
        "{file}.cram.crai"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule index_bam:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input}"


rule cram_to_bam:
    """
    A generic rule for converting cram to bam when downstream tools require bam.
    """
    input:
        cram="{file}.cram",
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
    output:
        temp("{file}.bam")
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools view -b {input.cram} -o {output} -T {input.fa}"

rule nsort_cram:
    """
    Generic rule for sorting a cram by read name.
    """
    input:
        crm="{file}.cram",
        fa=genome_fa,
        fai=genome_fai,
        gzi=genome_gzi,
    output:
        "{file}.nsrt.cram"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools sort -n -O cram {input.crm} --reference {input.fa} -o {output} "
        "--output-fmt-option lossy_names=1,level=9,store_md=0,store_nm=0"

# ------------------------------------------------------------------------------
# HELP
# ------------------------------------------------------------------------------

rule help:
    shell:
        """
        echo '
        ===== google cloud =====

        # install prereqs
        pip install kubernetes
        gcloud components install kubectl

        # set up cluster variables
        CLUSTER_NAME=snk-cl2
        NODES=8
        ZONE=us-central1-a
        REMOTE=GS
        PREFIX=archibald
        MACHINE_TYPE=n1-standard-2

        # initialize cluster
        gcloud container clusters create $CLUSTER_NAME \
            --num-nodes=$NODES \
            --scopes storage-rw \
            --machine-type=$MACHINE_TYPE \
            --zone $ZONE

        # register cluster info
        gcloud container clusters get-credentials $CLUSTER_NAME --zone $ZONE


        snakemake --kubernetes --use-conda \
            --default-remote-provider $REMOTE \
            --default-remote-prefix $PREFIX \
            --latency-wait 300 \
            --jobs 8 \
            --verbose \
            --debug-dag

        # shut down your cluster
        gcloud container clusters delete $CLUSTER_NAME --zone $ZONE

        '
        """

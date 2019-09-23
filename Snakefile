from os.path import realpath
from os.path import split as pathsplit
import subprocess
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


# info
__author__ = "Matt Lawlor"

# SETUP
shell.executable("/bin/bash")
GS = GSRemoteProvider()
FTP = FTPRemoteProvider()
configfile: "config.yaml"

# CONTSTANTS TODO integrate these
BT2_TRIM_SIZE = 30
MACS2_SHIFT_SIZE = 15


# REMOTE RESOURCES TODO deal with this elegantly - allow mouse as well
hg38_idx_files = ["seq-resources/NCBI-hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index." + x for x in ["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]]
hg38_idx_paths = [GS.remote(y) for y in hg38_idx_files]
hg38_fa = GS.remote("seq-resources/NCBI-hg38/hg38.fa.gz")
hg38_fai = GS.remote("seq-resources/NCBI-hg38/hg38.fa.gz.fai")
hg38_gzi = GS.remote("seq-resources/NCBI-hg38/hg38.fa.gz.gzi")
hg38_idx_pfx = "seq-resources/NCBI-hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
hg38_chroi_names = ["chr"+str(x) for x in range(1,23)] + ["chrX"]
hg38_bl = GS.remote("seq-resources/NCBI-hg38/ENCFF419RSJ.bed.gz")
hg38_gs = "2.7e9"


rule target:
    """
    Generate footprints, peaks, and a cut matrix from all samples in the manifest.
    """
    input:
        expand("{s}_peaks.mrg.bed",s=config.get("samples",None)),
        expand("{s}_fps.bed",s=config.get("samples",None))


# ------------------------------------------------------------------------------
# Preproc
# ------------------------------------------------------------------------------

#http://biolearnr.blogspot.com/2017/11/snakemake-using-inputoutput-values-in.html
rule align_bt2:
    """
    Align reads with bowtie2.
    """
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
        4
    params:
        idx_pfx = hg38_idx_pfx,
    shell:
        "bowtie2 --trim-to 3:30 --phred33 "
        "--no-discordant "
        "--no-unal "
        "-k 1 -X 800 "
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
        fa=hg38_fa,
        fai=hg38_fai,
        gzi=hg38_gzi,
        bl=hg38_bl,
        crai="{s}.raw.cram.crai"
    output:
        temp("{s}.bl.cram")
    conda:
        "environments/bedtools.yaml"
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
        fa=hg38_fa,
        fai=hg38_fai,
        gzi=hg38_gzi,
    output:
        temp("{s}.fixm.cram")
    conda:
        "environments/bowtie2.yaml"
    threads: 2
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
        fa=hg38_fa,
        fai=hg38_fai,
        gzi=hg38_gzi,
    output:
        "{s}.clean.cram"
    params:
        chr=hg38_chroi_names
    conda:
        "environments/bowtie2.yaml"
    threads:
        2
    shell:
        "samtools view -u -q 30 "
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
        #sum="{s}_summits.bed",
        #np="{s}_peaks.narrowPeak",
        mrg="{s}_peaks.mrg.bed"
    params:
        gs=hg38_gs
    conda:
        "environments/macs2.yaml"
    shadow: "shallow"
    threads:
        1
    shell:
        "macs2 callpeak -t {input} -f BAM " # TODO right now this only takes the 1st read in the pair...
        "--nomodel --shift -15 --keep-dup all "
        "-g {params.gs} -n {wildcards.s} --call-summits; "
        #"bedtools sort -i {wildcards.s}_peaks.narrowPeak | bedtools merge > {output.mrg}"
        #"mv {wildcards.s}_summits.bed {output.sum}; "
        #"mv  {output.np}"
        "bedtools sort -i {wildcards.s}_peaks.narrowPeak | bedtools merge > {output.mrg}; "

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
        "environments/pydnase.yaml"
    threads:
        4
    params:
        odir= "{s}-fp-tmp/"
    shell:
        "rm -rf {params.odir} && mkdir {params.odir} && "
        "wellington_footprints.py -p {threads} -o {wildcards.s} -A {input.bed} {input.bam} {params.odir}; "
        "mv {params.odir}/p\ value\ cutoffs/{wildcards.s}.WellingtonFootprints.-10.bed {output}; "
        "rm -rf {params.odir}"

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
        "environments/samtools.yaml"
    shell:
        "samtools index {input}"

rule index_bam:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    conda:
        "environments/samtools.yaml"
    shell:
        "samtools index {input}"


rule cram_to_bam:
    """
    A generic rule for converting cram to bam when downstream tools require bam.
    """
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

rule nsort_cram:
    """
    Generic rule for sorting a cram by read name.
    """
    input:
        crm="{file}.cram",
        fa=hg38_fa,
        fai=hg38_fai,
        gzi=hg38_gzi,
    output:
        "{file}.nsrt.cram"
    conda:
        "environments/samtools.yaml"
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
        NODES=2
        ZONE=us-central1-a
        REMOTE=GS
        PREFIX=archibald
        MACHINE_TYPE=n1-standard-4

        # initialize cluster
        gcloud container clusters create $CLUSTER_NAME
            --num-nodes=$NODES
            --scopes storage-rw
            --machine-type=$MACHINE_TYPE
            --zone $ZONE

        # register cluster info
        gcloud container clusters get-credentials $CLUSTER_NAME --zone $ZONE

        # run snakemake
        snakemake --kubernetes --use-conda
            --default-remote-provider $REMOTE
            --default-remote-prefix $PREFIX
            --latency-wait 300
            --jobs 2
            --keep-remote

        # shut down your cluster
        gcloud container clusters delete $CLUSTER_NAME --zone $ZONE

        '
        """

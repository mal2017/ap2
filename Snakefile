from os.path import realpath
from os.path import split as pathsplit
import subprocess

# info
__author__ = "Matt Lawlor"

# set shell
shell.executable("/bin/bash")


# rule to download from sra, or ena, or ftp
# rule to combine multiple fqs per end

rule target:
    input:
        "index.sr.mmi",
        expand("alignments/{s}.raw.sam",s=config.get("samples",None))


rule index_genome:
    input:
        config.get("genome",None)
    output:
        "index.sr.mmi"
    conda:
        "environments/minimap2.yaml"
    shell:
        "minimap2 -x sr -d {output} {input}"

rule align:
    input:
        r1 = lambda wc: config["samples"][wc.s]["fastq"]["r1"],
        r2 = lambda wc: config["samples"][wc.s]["fastq"]["r2"],
        idx = rules.index_genome.output,
    output:
        "alignments/{s}.raw.sam"
    conda:
        "environments/minimap2.yaml"
    threads:
        2
    shell:
        "minimap2 -ax sr {input.idx} -o {output} -t {threads} "
        "{input.r1} {input.r2}"

rule sam_to_bam:
    input:
        "{file}.sam"
    output:
        temp("{file}.bam")
    threads:
        2
    shell:
        "samtools sort {input} -@ {threads} -o {output}"

import os
import snakemake.io
import glob


FASTQ_DIR="/home/mhussien/pag/kSpider_Data"
OUTPUT_DIR="/home/mhussien/pag/kframes_scale_1000"
SAMPLES, = glob_wildcards(FASTQ_DIR + "/{sample}_1.fastq.gz")


rule all:
    input: expand(OUTPUT_DIR+"/{sample}.{ext}", sample=SAMPLES, ext=['phmap', 'extra'])

rule peFastQ_to_mqf:
    output:
        kf   = OUTPUT_DIR + "/{sample}.phmap",
        extra  = OUTPUT_DIR +  "/{sample}.extra",
    input:
        r1   = FASTQ_DIR + "/{sample}_1.fastq.gz",
        r2   = FASTQ_DIR + "/{sample}_2.fastq.gz",

    shell: "/home/mhussien/miniconda3/envs/pag/bin/kSpider kframes --r1 {input.r1} --r2 {input.r2} -k 25 --scale 1000"
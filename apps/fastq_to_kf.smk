import os
import snakemake.io
import glob

(SAMPLES,READS,) = glob_wildcards("raw/{sample}_{read}.fastq.gz")
READS=["1","2"]

FASQ_DIR = "/home/mhussien/pag/kSpider_Data/

rule all:
    input: expand("{FASQ_DIR}/{sample}.{ext}", sample=SAMPLES, ext=['mqf', 'extra'])

rule peFastQ_to_mqf:
    input:
        r1="{FASQ_DIR}/{sample}_1.fastq.gz",
        r2="{FASQ_DIR}/{sample}_2.fastq.gz"

    output: ["{sample}.mqf", "{sample}.extra"]

    # better to pass in the threads than to hardcode them in the shell command
    shell: "/home/mhussien/miniconda3/envs/pag/bin/kSpider kframes --r1 {input.r1} --r2 {input.r2} -k 25
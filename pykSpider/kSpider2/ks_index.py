#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import sys

import _kSpider_internal as kSpider_internal

import click
from kSpider2.click_context import cli
import os
import sys


@click.pass_context
def validate_names(ctx, names_file):
    '''validate names file for indexing'''

    ctx.obj.INFO("validating names file..")

    with open(names_file) as names:
        for i, line in enumerate(names, 1):
            if len(line.strip().split("\t")) != 2:
                ctx.obj.ERROR(
                    f"invalid names line detected at L{i}: '{line.strip()}'")


@cli.command(name="index_kmers", help_priority=1)
@click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
@click.option('-n', '--names', "names_file", required=True, type=click.Path(exists=True), help="Names file")
@click.option('-k', '--kmer-size', "kSize", required=True, type=click.IntRange(7, 31, clamp=False), help="kmer size")
@click.option('-c', '--chunk-size', "chunk_size", required=False, type=click.INT, default=3000, help="chunk size")
@click.option('--strand-specific', "strand_specific", is_flag=True)
@click.option('-o', '--output', "output_prefix", required=False, default=None, help="index output file prefix")
@click.pass_context
def kmers(ctx, fasta_file, names_file, kSize, output_prefix, chunk_size, strand_specific):
    '''FASTA file indexing by Kmers'''

    validate_names(names_file)

    if not output_prefix:
        output_prefix = os.path.basename(fasta_file)
        output_prefix = os.path.splitext(output_prefix)[0]
        output_prefix = "idx" + "_" + output_prefix

    ctx.obj.INFO("Indexing has begun, please wait ....")

    if strand_specific:
        kSpider_internal.index_kmers_nonCanonical(
            kSize, fasta_file, names_file, chunk_size, output_prefix)
    else:
        kSpider_internal.index_kmers(
            kSize, fasta_file, names_file, chunk_size, output_prefix)

    ctx.obj.SUCCESS("Indexing has completed.")


@cli.command(name="index_skipmers", help_priority=2)
@click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
@click.option('-n', '--names', "names_file", required=True, type=click.Path(exists=True), help="Names file")
@click.option('-k', '--kmer-size', "skipmers_kSize", required=True, type=click.INT, help="kmer size")
@click.option('-m', '--cycle-bases', "skipmers_m", required=True, type=click.INT, help="used bases per cycle")
@click.option('-n', '--cycle-length', "skipmers_n", required=True, type=click.INT, help="kmer size(cycle length")
@click.option('-c', '--chunk-size', "chunk_size", required=False, type=click.INT, default=3000, help="chunk size")
@click.option('-o', '--output', "output_prefix", required=False, default=None, help="index output file prefix")
@click.pass_context
def skipmers(ctx, fasta_file, names_file, skipmers_kSize, skipmers_m, skipmers_n, chunk_size, output_prefix):
    '''FASTA file indexing by Skipmers'''

    validate_names(names_file)

    if not output_prefix:
        output_prefix = os.path.basename(fasta_file)
        output_prefix = os.path.splitext(output_prefix)[0]
        output_prefix = "idx" + "_" + output_prefix

    if skipmers_n < 1 or skipmers_n < skipmers_m or skipmers_kSize < skipmers_m or skipmers_kSize % skipmers_m != 0:
        raise click.BadParameter(
            "Invalid skip-mer shape!\nConditions: 0 < m <= n < k & k must be multiple of m")

    ctx.obj.INFO("Indexing has begun, please wait ....")
    kSpider_internal.index_skipmers(
        skipmers_m, skipmers_n, skipmers_kSize, fasta_file, names_file, chunk_size, output_prefix)

    ctx.obj.SUCCESS("Indexing has completed.")


@cli.command(name="index_protein", help_priority=3)
@click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
@click.option('-n', '--names', "names_file", required=True, type=click.Path(exists=True), help="Names file")
@click.option('-k', '--kmer-size', "kSize", required=True, type=click.IntRange(7, 31, clamp=False), help="kmer size")
@click.option('-c', '--chunk-size', "chunk_size", required=False, type=click.INT, default=3000, help="chunk size")
@click.option('--dayhoff', "dayhoff", is_flag=True, show_default=True, default=False, help="use Dayhoff encoding")
@click.option('-o', '--output', "output_prefix", required=False, default=None, help="index output file prefix")
@click.pass_context
def protein(ctx, fasta_file, names_file, kSize, output_prefix, chunk_size, dayhoff):
    '''FASTA file indexing by Protein'''

    validate_names(names_file)

    if not output_prefix:
        output_prefix = os.path.basename(fasta_file)
        output_prefix = os.path.splitext(output_prefix)[0]
        output_prefix = "idx" + "_" + output_prefix

    ctx.obj.INFO("Indexing has begun, please wait ....")

    if dayhoff:
        kSpider_internal.index_dayhoff(
            kSize, fasta_file, names_file, chunk_size, output_prefix)
    else:
        kSpider_internal.index_dayhoff(
            kSize, fasta_file, names_file, chunk_size, output_prefix)

    ctx.obj.SUCCESS("Indexing has completed.")

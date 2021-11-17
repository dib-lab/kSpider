#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os


@cli.command(name="kframes", help_priority=6)
@click.option('-c', '--chunk-size', "chunk_size", required=False, type=click.INT, default=3000, help="chunk size")
@click.option('-k', '--kmer-size', "kSize", required=True, type=click.IntRange(7, 31, clamp=False), help="kmer size")
@click.option('-l', '--fastx-list', "fastx_list", type=click.Path(exists=True))
@click.option('--paired-end', "paired_end", is_flag=True, show_default=True, default=False, help="paired-end reads mode")
@click.option('--protein', "protein", is_flag=True, show_default=True, default=False, help="parsing protein")
@click.option('--dayhoff', "dayhoff", is_flag=True, show_default=True, default=False, help="parsing protein in dayhoff encoding")


@click.pass_context
def main(ctx, fastx_list, paired_end, chunk_size, kSize, protein, dayhoff):
    """
    Convert all files in a directory to kDataFrames
    """
    
    if protein and paired_end:
        ctx.obj.ERROR("Protein can't be paired-end.")

    if not paired_end:
        with open(fastx_list) as fastx:
            for line in fastx:
                file_path = line.strip()
                if not os.path.exists(file_path):
                    ctx.obj.ERROR(f"{file_path} file does not exist")
                
                ctx.obj.INFO(f"Processing {file_path}")
                if dayhoff:
                    kSpider_internal.protein_to_kDataFrame(file_path, kSize, chunk_size, True, os.path.basename(file_path))
                elif protein:
                    kSpider_internal.protein_to_kDataFrame(file_path, kSize, chunk_size, False, os.path.basename(file_path))
                else:
                    kSpider_internal.single_end_to_kDataFrame(file_path, kSize, chunk_size)

    ctx.obj.SUCCESS("All files has been converted.")

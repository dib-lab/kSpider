#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

from click.decorators import option
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os


@cli.command(name="kframes", help_priority=6)
@click.option('-c', '--chunk-size', "chunk_size", required=False, type=click.INT, default=3000, help="chunk size")
@click.option('-k', '--kmer-size', "kSize", required=True, type=click.IntRange(7, 31, clamp=False), help="kmer size")
@click.option('--fastx', "fastx", type=click.Path(exists=True), help = "FASTX file path, works with interleaved paired-end and protein", required= False)
@click.option('--r1', "r1", type=click.Path(exists=True), help = "paired-end FASTX R1 file", required= False)
@click.option('--r2', "r2", type=click.Path(exists=True), help = "paired-end FASTX R2 file", required= False)
@click.option('--protein', "protein", is_flag=True, show_default=True, default=False, help="parsing protein")
@click.option('--dayhoff', "dayhoff", is_flag=True, show_default=True, default=False, help="parsing protein in dayhoff encoding")
@click.pass_context
def main(ctx, fastx, r1, r2, chunk_size, kSize, protein, dayhoff):
    """
    Convert all files in a directory to kDataFrames
    """
    
    if protein and (r1 or r2):
        ctx.obj.ERROR("Protein can't be paired-end.")
    
    if fastx and (r1 or r2):
        ctx.obj.ERROR("You can use either --fastx or --r1 --r2.")
    
    if not fastx and not(r1 and r2):
        ctx.obj.ERROR("You need to provide --r1 --r2.")
   
    if protein and dayhoff:
        ctx.obj.ERROR("You can use either --protein or --dayhoff")
        
    paired_end_flag = False
    single_end_flag = False
    protein_flag = False
    dayhoff_flag = False
    
    if r1 or r2:
        paired_end_flag = True
    elif fastx and protein:
        protein_flag = True
    elif fastx and dayhoff:
        protein_flag = True
    elif fastx and not protein and not dayhoff:
        single_end_flag = True

    if protein_flag:
        kSpider_internal.protein_to_kDataFrame(fastx, kSize, chunk_size, False, os.path.basename(fastx))
    elif dayhoff_flag:
        kSpider_internal.protein_to_kDataFrame(fastx, kSize, chunk_size, True, os.path.basename(fastx))
    elif single_end_flag:
        kSpider_internal.single_end_to_kDataFrame(fastx, kSize, chunk_size)
    else:
        ctx.obj.ERROR("paired-end is not implemented yet")
                    
    ctx.obj.SUCCESS("File(s) has been converted.")
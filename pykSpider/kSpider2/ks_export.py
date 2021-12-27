#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import sys

import _kSpider_internal as kSpider_internal

import click
from kSpider2.click_context import cli

@cli.command(name="export", help_priority=5)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('--mapped-pairwise', "mapped_pairwise", is_flag=True, help="Write the pairwise TSV with original names", default = False)
@click.option('--dist-mat', "distance_matrix", is_flag=True, help="Convert pairwise matrix to nxn distance matrix", default = False)
@click.option('--newick', "newick", is_flag=True, help="Convert pairwise matrix to newick format", default = False)
@click.option('--percentage', "percentage", is_flag=True, help="Use %containment instead of number of kmers", default = False)
@click.pass_context
def main(ctx, index_prefix):
    """
    Convert kSpider pairwise to multiple formats.
    """
    ID_to_name = dict()    
    with open(index_prefix + ".namesMap", 'r') as namesMap:
        next(namesMap)  # skip the header
        for row in namesMap:
            row = row.strip().split()
            ID_to_name[int(row[0])] = row[1]


    ctx.obj.INFO("Constructing the containment pairwise matrix.")

    kSpider_internal.pairwise(index_prefix)

    ctx.obj.SUCCESS("Done.")

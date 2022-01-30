#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os
from glob import glob


@cli.command(name="index_datasets", help_priority=5)
@click.option('--dir', "sketches_dir", required = True, help="Sketches directory (must contain only the sketches)")
@click.pass_context
def main(ctx, sketches_dir):
    """
    Index all sketches in a directory.
    """
    if not os.path.exists(sketches_dir):
        ctx.obj.ERROR(f"{sketches_dir} does not exist!")
    
    all_extra = list(glob(f"{sketches_dir}/*extra"))
    all_sketches_phmap = glob(f"{sketches_dir}/*phmap")    
    all_sketches_mqf = glob(f"{sketches_dir}/*mqf")    
    
    if len(all_extra) != (len(all_sketches_phmap) + len(all_sketches_mqf)):
        ctx.obj.ERROR(f"Inconsistent sketches files.")
    
    ctx.obj.INFO(f"Indexing sketches in {sketches_dir}.")
    kSpider_internal.index_datasets(sketches_dir)

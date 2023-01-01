#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os
from glob import glob


@cli.command(name="index", help_priority=2)
@click.option('--dir', "sketches_dir", required=True, help="Sketches directory (must contain only the sketches)")
@click.option('-k', '--kmer-size', "kSize", required=False, default=0, type=click.INT, help="kmer size (only if using --sourmash)")
@click.option('--sourmash', "sourmash", is_flag=True, show_default=True, default=False, help="use sourmash sigs instead of kProcessor")
@click.pass_context
def main(ctx, sketches_dir, sourmash, kSize):
    """
    Index all sketches in a directory.
    """
    if not os.path.exists(sketches_dir):
        ctx.obj.ERROR(f"{sketches_dir} does not exist!")

    if sourmash:
        if not kSize:
            ctx.obj.ERROR(f"must select kSize when using --sourmash")
        ctx.obj.INFO(
            f"Indexing sourmash sigs in {sketches_dir} with kSize={kSize}.")
        kSpider_internal.sourmash_sigs_indexing(sketches_dir, kSize)
        ctx.obj.SUCCESS("DONE!")

    else:
        all_extra = list(glob(f"{sketches_dir}/*extra"))
        all_sketches_phmap = glob(f"{sketches_dir}/*phmap")
        all_sketches_mqf = glob(f"{sketches_dir}/*mqf")

        if len(all_extra) != (len(all_sketches_phmap) + len(all_sketches_mqf)):
            ctx.obj.ERROR(f"Inconsistent sketches files.")

        ctx.obj.INFO(f"Indexing sketches in {sketches_dir}.")
        kSpider_internal.index_datasets(sketches_dir)
        ctx.obj.SUCCESS("DONE!")

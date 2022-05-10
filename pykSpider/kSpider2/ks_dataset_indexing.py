#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os
from glob import glob


@cli.command(name="index_datasets", help_priority=5)
@click.option('--dir', "sketches_dir", required=True, help="Sketches directory (must contain only the sketches)")
@click.option('-k', '--kmer-size', "kSize", required=False, default=0, type=click.INT, help="kmer size (only if using --sourmash)")
@click.option('-s', '--scale', "scale", required=False, default=1, type=click.INT, help="scale (only if using --sourmash)")
@click.option('--sourmash', "sourmash", is_flag=True, show_default=True, default=False, help="use sourmash sigs instead of kProcessor")
@click.option('--fast', "fast", is_flag=True, show_default=True, default=False, help="if you're certain kSize is found in all smash sigs")
@click.pass_context
def main(ctx, sketches_dir, sourmash, kSize, fast, scale):
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
        if fast:
            kSpider_internal.sourmash_index_kp1_fast(sketches_dir, kSize, scale)
        else:
            kSpider_internal.sourmash_index_kp1(sketches_dir, kSize)

    else:
        all_extra = list(glob(f"{sketches_dir}/*extra"))
        all_sketches_phmap = glob(f"{sketches_dir}/*phmap")
        all_sketches_mqf = glob(f"{sketches_dir}/*mqf")

        if len(all_extra) != (len(all_sketches_phmap) + len(all_sketches_mqf)):
            ctx.obj.ERROR(f"Inconsistent sketches files.")

        ctx.obj.INFO(f"Indexing sketches in {sketches_dir}.")
        kSpider_internal.index_datasets(sketches_dir)

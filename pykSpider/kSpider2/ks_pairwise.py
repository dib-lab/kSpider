#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli


@cli.command(name="pairwise", help_priority=3)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-t', '--threads', "user_threads", default=1, required=False, type=int, help="number of cores")
@click.pass_context
def main(ctx, index_prefix, user_threads):
    """
    Generate containment pairwise matrix.
    """

    ctx.obj.INFO(
        f"Constructing the containment pairwise matrix using {user_threads} cores.")

    kSpider_internal.pairwise(index_prefix, user_threads)

    ctx.obj.SUCCESS("Done.")

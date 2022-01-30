#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import sys

import _kSpider_internal as kSpider_internal

import click
from kSpider2.click_context import cli


@cli.command(name="pairwise", help_priority=6)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.pass_context
def main(ctx, index_prefix):
    """
    Generate containment pairwise matrix.
    """

    ctx.obj.INFO("Constructing the containment pairwise matrix.")

    kSpider_internal.pairwise(index_prefix)

    ctx.obj.SUCCESS("Done.")

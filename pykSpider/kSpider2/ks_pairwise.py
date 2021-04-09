#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import sys

import _kSpider_internal as kSpider_internal

# try:
#     import pykSpider.kSpider_internal
# except ImportError:
#     print("kSpider_internal is not built yet", file = sys.stderr)

import click
from kSpider2.click_context import cli


@cli.command(name = "pairwise", help_priority=1)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="kProcessor index file prefix")
@click.pass_context
def main(ctx, index_prefix):
    """
    Generating containment pairwise matrix for kProcessor index.
    """

    kSpider_internal.pairwise(index_prefix)

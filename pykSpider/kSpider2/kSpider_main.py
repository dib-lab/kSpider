import sys
import click

from kSpider2.click_context import cli
from kSpider2.ks_pairwise import main as pairwise_main   # pylint: disable=relative-beyond-top-level

cli.add_command(pairwise_main, name="pairwise")

if __name__ == '__main__':
    cli()

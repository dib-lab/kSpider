import sys
import click

from kSpider.click_context import cli
from kSpider.ks_pairwise import main as pairwise_main   # pylint: disable=relative-beyond-top-level

cli.add_command(pairwise_main, name="pairwise")

if __name__ == '__main__':
    cli()

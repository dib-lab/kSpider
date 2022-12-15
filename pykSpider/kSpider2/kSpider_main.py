import sys

from kSpider2.click_context import cli
from kSpider2.ks_pairwise import main as pairwise_main   # pylint: disable=relative-beyond-top-level
# from kSpider2.ks_index import kmers, skipmers, protein   # pylint: disable=relative-beyond-top-level
from kSpider2.ks_clustering import main as clustering
from kSpider2.ks_fastx_to_kfs import main as fastx_to_kfs
from kSpider2.ks_export import main as export
from kSpider2.ks_dataset_indexing import main as index_datasets


# cli.add_command(kmers, name="index_kmers")
# cli.add_command(skipmers, name="index_skipmers")
# cli.add_command(protein, name="index_protein")
cli.add_command(fastx_to_kfs, name="sketch")
cli.add_command(index_datasets, name="index")
cli.add_command(pairwise_main, name="pairwise")
cli.add_command(clustering, name="cluster")
cli.add_command(export, name="export")


if __name__ == '__main__':
    cli()

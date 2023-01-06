from __future__ import division
from collections import defaultdict
import os
import click
from kSpider2.click_context import cli
import rustworkx as rx
from tqdm import tqdm


class Clusters:

    distance_to_col = {
        "min_cont": 3,
        "avg_cont": 4,
        "max_cont": 5,
        "ani": 6,
    }

    seq_to_kmers = dict()
    names_map = dict()

    def __init__(self, logger_obj, index_prefix, cut_off_threshold, dist_type):
        self.index_prefix = index_prefix
        self.Logger = logger_obj
        self.dist_type = dist_type
        self.edges_batch_number = 10_000_000
        self.names_file = index_prefix + ".namesMap"
        self.cut_off_threshold = cut_off_threshold
        self.seqToKmers_file = index_prefix + "_kSpider_seqToKmersNo.tsv"
        self.pairwise_file = index_prefix + "_kSpider_pairwise.tsv"
        self.output = index_prefix + \
            f"_kSpider_clusters_{cut_off_threshold}%.tsv"
        self.shared_kmers_threshold = 200
        self.Logger.INFO("Loading TSV pairwise file")
        self.load_seq_to_kmers(self.seqToKmers_file)
        self.tsv_get_namesmap()
        if dist_type not in self.distance_to_col:
            logger_obj.ERROR("unknown distance!")
        self.dist_col = self.distance_to_col[dist_type]

        if dist_type == "ani":
            if not os.path.exists(self.index_prefix + "_kSpider_pairwise.ani_col.tsv"):
                logger_obj.ERROR(f"ANI was selected, but the ani file {self.index_prefix}_kSpider_pairwise.ani_col.tsv was not found!")

        self.graph = rx.PyGraph()
        self.nodes_indeces = self.graph.add_nodes_from(
            list(self.names_map.keys()))

    def load_seq_to_kmers(self, tsv):
        with open(tsv) as KMER_COUNT:
            next(KMER_COUNT)
            for line in KMER_COUNT:
                seq_ID, no_of_kmers = tuple(line.strip().split('\t')[1:])
                self.seq_to_kmers[int(seq_ID)] = int(no_of_kmers)

    def tsv_get_namesmap(self):
        with open(self.names_file, 'r') as namesMap:
            next(namesMap)  # skip the header
            for row in namesMap:
                row = row.strip().split()
                self.names_map[int(row[0])] = row[1]

    def construct_graph(self):
        batch_counter = 0
        edges_tuples = []

        print("[i] constructing graph")
        with open(self.pairwise_file, 'r') as pairwise_tsv:
            next(pairwise_tsv)  # skip header
            if self.dist_type == "ani":
                with open(self.index_prefix + "_kSpider_pairwise.ani_col.tsv") as ani_col_file:
                    for row in pairwise_tsv:
                        row = row.strip().split('\t')
                        seq1 = int(row[0]) - 1
                        seq2 = int(row[1]) - 1
                        distance = float(next(ani_col_file).strip()) * 100

                        # don't make graph edge
                        if distance < self.cut_off_threshold:
                            continue

                        if batch_counter < self.edges_batch_number:
                            batch_counter += 1
                            edges_tuples.append((seq1, seq2, distance))
                        else:
                            self.graph.add_edges_from(edges_tuples)
                            batch_counter = 0
                            edges_tuples.clear()

                    else:
                        if len(edges_tuples):
                            self.graph.add_edges_from(edges_tuples)
            else:
                for row in pairwise_tsv:
                    row = row.strip().split('\t')
                    seq1 = int(row[0]) - 1
                    seq2 = int(row[1]) - 1
                    distance = float(row[self.dist_col]) * 100

                    # don't make graph edge
                    if distance < self.cut_off_threshold:
                        continue

                    if batch_counter < self.edges_batch_number:
                        batch_counter += 1
                        edges_tuples.append((seq1, seq2, distance))
                    else:
                        self.graph.add_edges_from(edges_tuples)
                        batch_counter = 0
                        edges_tuples.clear()

                else:
                    if len(edges_tuples):
                        self.graph.add_edges_from(edges_tuples)

    def cluster_graph(self):
        self.connected_components = rx.connected_components(self.graph)
        self.Logger.INFO(
            f"number of clusters: {len(self.connected_components)}")
        single_components = 0
        retworkx_export = self.index_prefix + \
            f"_kSpider_graph_{self.cut_off_threshold}%.json"
        # and {self.output} ...")
        self.Logger.INFO(f"writing {retworkx_export}")
        # rx.node_link_json(self.graph, path = retworkx_export)
        with open(self.output, 'w') as CLUSTERS:
            for component in self.connected_components:
                # uncomment to exclude single genome clusters from exporting
                # if len(component) == 1:
                #     single_components += 1
                #     continue
                named_component = [self.names_map[node + 1]
                                   for node in component]
                CLUSTERS.write(','.join(named_component) + '\n')


"""
TODO:
New help messages

1. containment cutoff (sim_cutoff): cluster sequences with (containment > cutoff) where containment = shared kmers % to the total kmers in the smallest node.
2. connectivity cutoff (con_cutoff): cluster sequences with (connectivity > cutoff) where connectivity = shared kmers % to the total kmers in the largest node.
3. min count cutoff (min_count): the min kmers count of a node to connect two clusters, otherwise the node will be reported twice in both clusters.
"""


@cli.command(name="cluster", help_priority=4)
@click.option('-c', '--cutoff', required=False, type=click.FloatRange(0, 1, clamp=False), default=0.0, show_default=True, help="cluster sequences with (containment > cutoff)")
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING, help="Index file prefix")
@click.option('-d', '--dist-type', "distance_type", required=False, default="max_cont", show_default=True, type=click.STRING, help="select from ['min_containment', 'avg_containment', 'max_containment', 'ani']")
@click.pass_context
def main(ctx, index_prefix, cutoff, distance_type):
    """Sequence clustering."""
    cutoff = float(cutoff) * 100
    kCl = Clusters(logger_obj=ctx.obj, index_prefix=index_prefix,
                   cut_off_threshold=cutoff, dist_type=distance_type)
    ctx.obj.INFO("Building the main graph...")
    kCl.construct_graph()
    ctx.obj.INFO("Clustering...")
    kCl.cluster_graph()

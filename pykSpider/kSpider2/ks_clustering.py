from __future__ import division
from collections import defaultdict
import itertools
import sys
import os
import sqlite3
import click
from kSpider2.click_context import cli
import glob


class kClusters:

    source = []
    target = []
    source2 = []
    target2 = []
    seq_to_kmers = dict()
    names_map = dict()
    components = defaultdict(set)

    def __init__(self, logger_obj, index_prefix, cut_off_threshold):
        self.Logger = logger_obj
        self.names_file = index_prefix + ".namesMap"
        self.cut_off_threshold = cut_off_threshold
        self.seqToKmers_file = index_prefix + "_kSpider_seqToKmersNo.tsv"
        self.pairwise_file = index_prefix + "_kSpider_pairwise.tsv"
        self.uncovered_seqs = set()
        self.shared_kmers_threshold = 200
        self.seq_to_clusterid = dict()
        self.max_cluster_id = 0
        self.Logger.INFO("Loading TSV pairwise file")
        self.load_seq_to_kmers(self.seqToKmers_file)
        self.tsv_get_namesmap()

    def load_seq_to_kmers(self, tsv):
        with open(tsv) as KMER_COUNT:
            next(KMER_COUNT)
            for line in KMER_COUNT:
                seq_ID, no_of_kmers = tuple(line.strip().split('\t')[1:])
                self.seq_to_kmers[int(seq_ID)] = int(no_of_kmers)

    def ids_to_names(self, cluster):
        new_cluster = []
        for _id in cluster:
            new_cluster.append(self.names_map[int(_id)])

        return new_cluster

    def tsv_get_namesmap(self):
        with open(self.names_file, 'r') as namesMap:
            next(namesMap)  # skip the header
            for row in namesMap:
                row = row.strip().split()
                self.names_map[int(row[0])] = row[1]

    def tsv_build_graph(self):

        with open(self.pairwise_file, 'r') as pairwise_tsv:
            next(pairwise_tsv)  # skip header

            for row in pairwise_tsv:
                row = row.strip().split()
                seq1 = int(row[1])
                seq2 = int(row[2])
                shared_kmers = int(row[3])
                containment = 0.0

                min_seq = float(
                    min(self.seq_to_kmers[seq1], self.seq_to_kmers[seq2]))
                containment = shared_kmers / min_seq

                if containment < self.cut_off_threshold:
                    continue

                if shared_kmers < self.shared_kmers_threshold:
                    self.source2.append(seq1)
                    self.target2.append(seq2)

                elif shared_kmers >= self.shared_kmers_threshold:
                    self.source.append(seq1)
                    self.target.append(seq2)

            # # For covering clusters with single sequence
            uncovered_seqs_1 = set(self.names_map.keys()) - \
                set(self.source).union(set(self.target))
            for seq in uncovered_seqs_1:
                self.uncovered_seqs.add(seq)

            # OR:
            # for i in range(1, len(self.names_map) + 1, 1):
            #     self.source.append(i)
            #     self.target.append(i)

    def clustering(self):
        registers = defaultdict(lambda: None)

        def find(x):
            l = registers[x]
            if l is not None:
                l = find(l)
                registers[x] = l
                return l
            return x

        def union(x, y):
            lx, ly = find(x), find(y)
            if lx != ly:
                registers[lx] = ly

        for i in range(len(self.source)):
            union(self.source.pop(), self.target.pop())

        for x in registers:
            self.components[find(x)].add(x)

        temp_components = self.components.copy()
        self.components.clear()

        for cluster_id, (k, v) in enumerate(temp_components.items(), 1):
            self.components[cluster_id] = set(v)
            for seq in v:
                self.seq_to_clusterid[seq] = cluster_id

        temp_components.clear()
        self.post_clustering()

    def post_clustering(self):
        registers2 = defaultdict(lambda: None)
        local_components = defaultdict(set)
        covered_seqs = set()

        def find(x):
            l = registers2[x]
            if l is not None:
                l = find(l)
                registers2[x] = l
                return l
            return x

        def union(x, y):
            lx, ly = find(x), find(y)
            if lx != ly:
                registers2[lx] = ly

        for i in range(len(self.source2)):
            union(self.source2.pop(), self.target2.pop())

        for x in registers2:
            local_components[find(x)].add(x)

        self.components = dict(self.components)

        covered_clusters = set()

        for cluster2_id, (k, v) in enumerate(local_components.items(), 1):

            for seq in v:
                covered_seqs.add(seq)

            for seq in v:
                if seq in self.seq_to_clusterid:
                    cluster_id = self.seq_to_clusterid[seq]
                    to_be_added = set()

                    for i in v:
                        if i not in self.seq_to_clusterid:
                            to_be_added.add(i)

                    self.components[cluster_id] = self.components[cluster_id].union(
                        to_be_added)
                    covered_clusters.add(k)
                    continue

        self.uncovered_seqs = self.uncovered_seqs - covered_seqs
        uncovered_clusters = set(local_components.keys()) - covered_clusters
        max_id = len(self.components)
        for i, unc in enumerate(uncovered_clusters, 1):
            max_id += 1
            self.components[max_id] = local_components[unc]

        for seq in self.uncovered_seqs:
            max_id += 1
            self.components[max_id] = {seq}

    def export_kCluster(self):
        kCluster_file_name = f"kSpider_{self.cut_off_threshold:.2f}%_"
        kCluster_file_name += os.path.basename(
            self.pairwise_file).split(".")[0]
        kCluster_file_name += ".clusters.tsv"

        with open(kCluster_file_name, 'w') as kClusters:
            kClusters.write("kClust_id\tseqs_ids\n")
            for cluster_id, (k, v) in enumerate(self.components.items(), 1):
                kClusters.write(
                    f"{cluster_id}\t{'|'.join(self.ids_to_names(v))}\n")

        self.Logger.INFO(f"Total Number Of Clusters: {cluster_id}")


"""
TODO:
New help messages

1. containment cutoff (sim_cutoff): cluster sequences with (containment > cutoff) where containment = shared kmers % to the total kmers in the smallest node.
2. connectivity cutoff (con_cutoff): cluster sequences with (connectivity > cutoff) where connectivity = shared kmers % to the total kmers in the largest node.
3. min count cutoff (min_count): the min kmers count of a node to connect two clusters, otherwise the node will be reported twice in both clusters.
"""


@cli.command(name="cluster", help_priority=7)
@click.option('-c', '--cutoff', required=False, type=click.FloatRange(0, 1, clamp=False), default=0.0, show_default=True, help="cluster sequences with (containment > cutoff)")
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING, help="Index file prefix")
@click.pass_context
def main(ctx, index_prefix, cutoff):
    """Sequence clustering."""

    kCl = kClusters(logger_obj=ctx.obj,
                    index_prefix=index_prefix, cut_off_threshold=cutoff)
    ctx.obj.INFO("Building the main graph...")
    kCl.tsv_build_graph()
    ctx.obj.INFO("Clustering...")
    kCl.clustering()
    ctx.obj.INFO("Exporting ...")
    kCl.export_kCluster()

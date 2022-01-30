#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import os
import _kSpider_internal as kSpider_internal
import click
import pandas as pd
from kSpider2.click_context import cli
from scipy.cluster.hierarchy import linkage, to_tree, ClusterWarning
from warnings import simplefilter
simplefilter("ignore", ClusterWarning)


# Thanks to https://stackoverflow.com/a/31878514/3371177
def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist,
                            leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist,
                            leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick


@cli.command(name="export", help_priority=8)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('--dist-mat', "distance_matrix", is_flag=True, help="Convert pairwise matrix to nxn distance matrix", default=False)
@click.option('--newick', "newick", is_flag=True, help="Convert pairwise (containment) matrix to newick format", default=False)
@click.option('--containment', "containment", is_flag=True, help="Use max %containment instead of number of kmers", default=False)
@click.pass_context
def main(ctx, index_prefix, containment, newick, distance_matrix):
    """
    Export kSpider pairwise to multiple formats.
    """

    index_basename = os.path.basename(index_prefix)
    kSpider_pairwise_tsv = f"{index_prefix}_kSpider_pairwise.tsv"
    namesMap_file = f"{index_prefix}.namesMap"
    seqToKmers_tsv = f"{index_prefix}_kSpider_seqToKmersNo.tsv"
    # Check for existing pairwise file
    for _file in [kSpider_pairwise_tsv, namesMap_file, seqToKmers_tsv]:
        if not os.path.exists(_file):
            ctx.obj.ERROR("File {_file} is not found.")

    """
    # Load kmer count per record
    """
    seq_to_kmers = dict()
    with open(seqToKmers_tsv) as KMER_COUNT:
        next(KMER_COUNT)
        for line in KMER_COUNT:
            seq_ID, no_of_kmers = tuple(line.strip().split('\t')[1:])
            seq_to_kmers[seq_ID] = int(no_of_kmers)

    """
    # Parse namesmap
    """

    namesMap_dict = dict()
    with open(namesMap_file) as NAMES:
        next(NAMES)
        for line in NAMES:
            line = line.strip().split()
            _id = line[0]
            _name = line[1]
            namesMap_dict[_id] = _name

    """Parse kSpider's pairwise
    """
    distances = dict()
    labeled_out = f"kSpider_{index_basename}_pairwise.tsv"
    distmatrix_out = f"kSpider_{index_basename}_distmat.tsv"
    newick_out = f"kSpider_{index_basename}.newick"

    with open(kSpider_pairwise_tsv) as PAIRWISE, open(labeled_out, 'w') as NEW:
        # Skip header
        third_column = "shared_kmers" if not containment else "max_containment"
        ctx.obj.INFO(f"Writing pairwise matrix to {labeled_out}")
        NEW.write(f"grp1\tgrp2\t{third_column}\n")
        next(PAIRWISE)
        for line in PAIRWISE:
            line = (line.strip().split('\t')[1:])
            origin_grp1 = line[0]
            origin_grp2 = line[1]
            shared_kmers = int(line[2])
            grp1 = namesMap_dict[origin_grp1]
            grp2 = namesMap_dict[origin_grp2]
            min_seq = float(
                min(seq_to_kmers[origin_grp1], seq_to_kmers[origin_grp2]))
            max_containment = (shared_kmers / min_seq)
            value = shared_kmers
            if containment:
                value = max_containment
            if distance_matrix or newick:
                distances[(grp1, grp2)] = max_containment

            NEW.write(f"{grp1}\t{grp2}\t{value}\n")

    if distance_matrix or newick:
        unique_ids = sorted(set([x for y in distances.keys() for x in y]))
        df = pd.DataFrame(index=unique_ids, columns=unique_ids)
        for k, v in distances.items():
            df.loc[k[0], k[1]] = 1-v
            df.loc[k[1], k[0]] = 1-v

        df = df.fillna(0)
        if distance_matrix:
            ctx.obj.INFO(f"Writing distance matrix to {distmatrix_out}")
            df.to_csv(distmatrix_out, sep='\t')
        if newick:
            if not distance_matrix:
                df.to_csv(distmatrix_out, sep='\t')
            loaded_df = pd.read_csv(distmatrix_out, sep='\t')
            os.remove(distmatrix_out)
            ctx.obj.INFO(f"Writing newick to {newick_out}.")
            names = list(loaded_df.columns[1:])
            dist = loaded_df[loaded_df.columns[1:]].to_numpy()
            Z = linkage(dist, 'complete')
            tree = to_tree(Z, False)

            newick = get_newick(tree, tree.dist, names)
            with open(newick_out, 'w') as NW:
                NW.write(newick)

    ctx.obj.SUCCESS("Done.")

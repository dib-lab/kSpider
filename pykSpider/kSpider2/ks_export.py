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


@cli.command(name="export", help_priority=5)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
# @click.option('--dist-mat', "distance_matrix", is_flag=True, help="Convert pairwise matrix to NxN distance matrix", default=False)
@click.option('--newick', "newick", is_flag=True, help="Convert pairwise (containment) matrix to newick format", default=False)
@click.option('-d', '--dist-type', "distance_type", required=False, default="max_cont", show_default=True, type=click.STRING, help="select from ['min_containment', 'avg_containment', 'max_containment', 'ani']")
@click.option('-o', "overwritten_output", default="na", required=False, type=click.STRING, help="custom output file name prefix")
@click.pass_context
def main(ctx, index_prefix, newick, distance_type, overwritten_output):
    """
    Export kSpider pairwise to multiple formats.
    """
    
    index_basename = os.path.basename(index_prefix)
    kSpider_pairwise_tsv = f"{index_prefix}_kSpider_pairwise.tsv"
    namesMap_file = f"{index_prefix}.namesMap"
    seqToKmers_tsv = f"{index_prefix}_kSpider_seqToKmersNo.tsv"
    
    LOGGER = ctx.obj
    
    distance_to_col = {
        "min_cont": 3,
        "avg_cont": 4,
        "max_cont": 5,
        "ani": 99
    }
    
    if distance_type not in distance_to_col:
        LOGGER.ERROR("unknown distance!")
    
    dist_col = distance_to_col[distance_type]
    if dist_col == "ani":
        with open(kSpider_pairwise_tsv, 'r') as pairwise_tsv:
            if "ani" not in next(pairwise_tsv).lower():
                LOGGER.ERROR("ANI was selected but was not found in the pairwise file.\nPlease, run kSpider pairwise --extend_with_ani -i <index_prefix> script")
    
    
    
    # Check for existing pairwise file
    for _file in [kSpider_pairwise_tsv, namesMap_file, seqToKmers_tsv]:
        if not os.path.exists(_file):
            LOGGER.ERROR(f"File {_file} is not found.")

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
    
    if overwritten_output != "na":
        labeled_out = f"{overwritten_output}_pairwise.tsv"
        distmatrix_out = f"{overwritten_output}_distmat.tsv"
        newick_out = f"{overwritten_output}.newick"
    
    if distance_type == "ani":
        with open(kSpider_pairwise_tsv) as PAIRWISE, open(labeled_out, 'w') as NEW, open(index_prefix + "_kSpider_pairwise.ani_col.tsv") as ANI:
            ctx.obj.INFO(f"Writing pairwise matrix to {labeled_out}")
            NEW.write(f"grp1\tgrp2\t{distance_type}\n")
            # Skip header
            next(PAIRWISE)
            next(ANI)
            for line in PAIRWISE:
                line = (line.strip().split('\t'))
                origin_grp1 = line[0]
                origin_grp2 = line[1]
                grp1 = namesMap_dict[origin_grp1]
                grp2 = namesMap_dict[origin_grp2]
                dist_metric = float(next(ANI).strip())
                distances[(grp1, grp2)] = dist_metric
                NEW.write(f"{grp1}\t{grp2}\t{dist_metric}\n")
                
    else:
        with open(kSpider_pairwise_tsv) as PAIRWISE, open(labeled_out, 'w') as NEW:
            ctx.obj.INFO(f"Writing pairwise matrix to {labeled_out}")
            NEW.write(f"grp1\tgrp2\t{distance_type}\n")
            # Skip header
            next(PAIRWISE)
            for line in PAIRWISE:
                line = (line.strip().split('\t'))
                origin_grp1 = line[0]
                origin_grp2 = line[1]
                grp1 = namesMap_dict[origin_grp1]
                grp2 = namesMap_dict[origin_grp2]
                dist_metric = float(line[dist_col])
                distances[(grp1, grp2)] = dist_metric
                NEW.write(f"{grp1}\t{grp2}\t{dist_metric}\n")

    unique_ids = sorted(set([x for y in distances.keys() for x in y]))
    df = pd.DataFrame(index=unique_ids, columns=unique_ids)
    for k, v in distances.items():
        df.loc[k[0], k[1]] = 1-v
        df.loc[k[1], k[0]] = 1-v

    df = df.fillna(0)
    LOGGER.INFO(f"Writing distance matrix to {distmatrix_out}")
    df.to_csv(distmatrix_out, sep='\t')
        
    if newick:
        loaded_df = pd.read_csv(distmatrix_out, sep='\t')
        LOGGER.INFO(f"Writing newick to {newick_out}.")
        names = list(loaded_df.columns[1:])
        dist = loaded_df[loaded_df.columns[1:]].to_numpy()
        Z = linkage(dist, 'single')
        tree = to_tree(Z, False)

        newick = get_newick(tree, tree.dist, names)
        with open(newick_out, 'w') as NW:
            NW.write(newick)

    LOGGER.SUCCESS("Done.")

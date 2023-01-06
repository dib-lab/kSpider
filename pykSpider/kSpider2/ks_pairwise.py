#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli



@cli.command(name="pairwise", help_priority=3)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('--estimate-ani', "ani", is_flag=True, show_default=True, default=False, help="estimate ANI and write result in a new file with single column")
@click.option('-t', '--threads', "user_threads", default=1, required=False, type=int, help="number of cores")
@click.option('-s', '--scale', "sourmash_scale", required=False, default=0, type=int, help="scale used in creating sourmash sigs (only when using --estimate-ani)")
@click.pass_context
def main(ctx, index_prefix, user_threads, ani, sourmash_scale):
    """
    Generate containment pairwise matrix.
    """
    if not ani:
        ctx.obj.INFO(
            f"Constructing the containment pairwise matrix using {user_threads} cores.")
        kSpider_internal.pairwise(index_prefix, user_threads)
        ctx.obj.SUCCESS("Done.")
    else:
        if user_threads > 1:
            ctx.obj.WARNING("sorry, current ANI estimation does not allow multithreading")
            
        if not sourmash_scale:
            ctx.obj.ERROR("estimating ANI requires to provide --scale value")
        
        # Detect kmer size
        kSize = None
        with open(f"{index_prefix}.extra") as extra:
            kSize = int(next(extra))
        
        from sourmash.distance_utils import containment_to_distance
        pairwise_file = index_prefix + "_kSpider_pairwise.tsv"
        ani_col = index_prefix + "_kSpider_pairwise.ani_col.tsv"
        seqToKmers_tsv = f"{index_prefix}_kSpider_seqToKmersNo.tsv"
        
        
        """
        # Load kmer count per record
        """
        id_to_kmer_count = dict()
        with open(seqToKmers_tsv) as KMER_COUNT:
            next(KMER_COUNT)
            for line in KMER_COUNT:
                seq_ID, no_of_kmers = tuple(line.strip().split('\t')[1:])
                id_to_kmer_count[int(seq_ID)] = int(no_of_kmers)
        
        
        
        with open(pairwise_file) as PAIRWISE, open(ani_col, 'w') as ANI_COL:
            next(PAIRWISE)
            ANI_COL.write("avg_ani\n")
            for or_line in PAIRWISE:
                line = or_line.strip().split('\t')
                shared_kmers = int(line[2])
                if shared_kmers < 5:
                    continue
                id_1 = int(line[0])
                id_2 = int(line[1])
                min_containment = float(line[3])
                max_containment = float(line[5])
                # print("Debug:")
                # print(
                #     (min_containment, kSize, sourmash_scale, id_to_kmer_count[]*sourmash_scale)
                # )
                ani_1_in_2 = containment_to_distance(min_containment, kSize, sourmash_scale, n_unique_kmers= id_to_kmer_count[id_2]*sourmash_scale).ani
                ani_2_in_1 = containment_to_distance(max_containment, kSize, sourmash_scale, n_unique_kmers= id_to_kmer_count[id_1]*sourmash_scale).ani
                avg_ani = (ani_1_in_2 + ani_2_in_1) / 2.0       
                new_line = f"{avg_ani}\n"
                ANI_COL.write(new_line)

        
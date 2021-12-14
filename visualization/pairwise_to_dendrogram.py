import pandas as pd
import sys
import os
import re
import plotly.figure_factory as ff
import numpy as np


if len(sys.argv) <5:
    sys.exit("run: python pairwise_to_containment_distance.py <pairwise_tsv> <seq_to_kmers_tsv> <namesMap_file> <sraTable_csv>")


pairwise_file = sys.argv[1]
seqToKmers_tsv = sys.argv[2]
namesMap_file = sys.argv[3]
sraTable_csv = sys.argv[4]

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
# Convert IDs to meaningful names
"""

run_ID_OUTPUT = dict()
with open(sraTable_csv) as RUN_TABLE:
    next(RUN_TABLE) # skip header
    for line in RUN_TABLE:
        line = line.strip()
        line = re.split(r',(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)', line)
        run_id = line[0]
        biosample = line[3]
        biome = line[16]
        geolocation = line[21]
        run_ID_OUTPUT[run_id] = f"{biosample}_{biome}_{geolocation}"


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



obj_distance = dict()
with open(pairwise_file) as PAIRWISE:
    # Skip header
    next(PAIRWISE)
    for line in PAIRWISE:
        line = (line.strip().split('\t')[1:])
        sample1_old, sample2_old, shared_kmers = line
        sample1 = namesMap_dict[sample1_old] + '_' + run_ID_OUTPUT[namesMap_dict[sample1_old]]
        sample2 = namesMap_dict[sample2_old] + '_' + run_ID_OUTPUT[namesMap_dict[sample2_old]]
        min_seq = float(min(seq_to_kmers[sample1_old], seq_to_kmers[sample2_old]))
        containment = int(shared_kmers) / min_seq
        obj_distance[(sample1, sample2)] = containment

unique_ids = sorted(set([x for y in obj_distance.keys() for x in y]))
df = pd.DataFrame(index=unique_ids, columns=unique_ids)

for k, v in obj_distance.items():
    df.loc[k[0], k[1]] = v
    df.loc[k[1], k[0]] = v

df = df.fillna(0)

df.to_csv("distmat_" + os.path.basename(pairwise_file), sep= '\t')

df = pd.read_csv("distmat_" + os.path.basename(pairwise_file), sep = '\t')
names = list(df.columns[1:])
distances = df[df.columns[1:]].to_numpy()
fig = ff.create_dendrogram(distances, color_threshold=0.1, labels = names)
fig.update_layout(width=4500, height=1000)
fig.write_html("wide_cont_dendro.html")

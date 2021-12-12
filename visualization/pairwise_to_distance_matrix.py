import pandas as pd
import sys
import os
import re

if len(sys.argv) < 2:
    sys.exit("run: python pairwise_to_distane_matrix.py <pairwise_tsv>")

pairwise_file = sys.argv[1]

obj_distance = dict()
with open(pairwise_file) as PAIRWISE:
    # Skip header
    next(PAIRWISE)
    for line in PAIRWISE:
        line = (line.strip().split('\t')[1:])
        sample1, sample2, shared_kmers = line
        obj_distance[(sample1, sample2)] = int(shared_kmers)

unique_ids = sorted(set([x for y in obj_distance.keys() for x in y]))
df = pd.DataFrame(index=unique_ids, columns=unique_ids)

for k, v in obj_distance.items():
    df.loc[k[0], k[1]] = v
    df.loc[k[1], k[0]] = v

df = df.fillna(0)

df.to_csv("distmat_" + os.path.basename(pairwise_file), sep= '\t')

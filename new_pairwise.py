import retworkx as rx
from tqdm import tqdm
from itertools import combinations as comb
from collections import defaultdict
import json

namesmap_file = "/home/mabuelanin/dib-dev/nasa/kSpider_optimize/work_dir/subset_sigs_100k.namesMap"
intvectors_file = "/home/mabuelanin/dib-dev/nasa/kSpider_optimize/work_dir/subset_sigs_100kcolors.intvectors"
colors_count_file = "/home/mabuelanin/dib-dev/nasa/kSpider_optimize/work_dir/subset_sigs_100k_kSpider_colorCount.tsv"

vertices = []
color_to_sources = dict()
color_to_count = dict()

source_to_kmer_count = dict()

# 1 load namesmap to construct vertices
"""
1753
1716 SRR9224419
1025 SRR5580024
"""
print("[i] parsing namesmap file")
with open(namesmap_file) as NAMES:
    total_sources = int(next(NAMES).strip())
    for line in tqdm(NAMES, total = total_sources):
        source_id = int(line.split(' ')[0])
        vertices.append(source_id)
        source_to_kmer_count[source_id] = 0 # related to #4

# 2 load colors and sources
"""
353893
217587 2 70 1132
162410 5 671 816 861 1078 1082
"""
print("[i] loading colors and sources")
with open(intvectors_file) as INTVECTORS:
    total_colors = int(next(INTVECTORS).strip())
    for line in tqdm(INTVECTORS, total=total_colors):
        line = list(map(int, line.strip().split(' ')))
        color_to_sources[line[0]] = line[2:]

# 3 load colors count
"""
color,count
217587,3
162410,1
"""
print("[i] loading colors and count")
with open(colors_count_file) as CCOUNT:
    next(CCOUNT)
    for line in tqdm(CCOUNT, total = len(color_to_sources)):
        line = list(map(int, line.strip().split(',')))
        color_to_count[line[0]] = line[1]

# 4 kmer counting
print("[i] counting kmers")
for color, sources in tqdm(color_to_sources.items()):
    if color not in color_to_count:
        color_to_count[color] = 0
        continue
    ccount = color_to_count[color]
    for source in sources:
        source_to_kmer_count[source] += ccount

# 5 pairwise construction
graph = rx.PyGraph()
nodes_indeces = graph.add_nodes_from(vertices)


def def_value(): return 0


batch_size = 50000000
batch_counter = 0
edges_batch = defaultdict(def_value)
edges_tuples = []

print("[i] constructing graph")
for color, sources in tqdm(color_to_sources.items()):
    for pair in comb(sources, 2):
        if batch_counter < batch_size:
            batch_counter += 1            
            p1 = pair[0] - 1
            p2 = pair[1] - 1
            edges_batch[(p1,p2)] += color_to_count[color]
        else:
            for (p1, p2), weight in edges_batch.items():
                try:
                    x = graph.get_edge_data(p1, p2)
                    edges_tuples.append((p1, p2, x + weight))
                except:
                    edges_tuples.append((p1, p2, weight))

            graph.add_edges_from(edges_tuples)
            batch_counter = 0
            edges_batch.clear()
            edges_tuples.clear()

# print all edges
all_edges = graph.edge_index_map()
with open("graph.json", "w") as outfile:
    json.dump(dict(all_edges), outfile, indent = 4)
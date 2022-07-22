import retworkx as rx
from tqdm import tqdm
import sys

if len(sys.argv) != 3:
    sys.exit("run: python rx_clustering.py <index_prefix> <clustering_threshold%>")

index_prefix = sys.argv[1]
CONTAINMENT_THRESHOLD = sys.argv[2]

namesmap_file = index_prefix + ".namesMap"
intvectors_file = index_prefix + "colors.intvectors"
colors_count_file = index_prefix + "_kSpider_colorCount.tsv"
pairwise_file = index_prefix + "_kSpider_pairwise.tsv"
kmer_count_file = index_prefix + "_kSpider_seqToKmersNo.tsv"

CONTAINMENT_THRESHOLD = 80

vertices = []
color_to_sources = dict()
color_to_count = dict()

source_to_kmer_count = dict()
id_to_name = dict()

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
        line = line.strip().split(' ')
        source_id = int(line[0])
        source_name = line[1]
        vertices.append(source_id)
        id_to_name[source_id] = source_name

# 2 kmer counting
print("[i] counting kmers")
with open(kmer_count_file) as KCOUNT_FILE:
    next(KCOUNT_FILE)
    for line in KCOUNT_FILE:
        line = line.strip().split()
        seq_id = int(line[1])
        number_of_kmers = int(line[2])
        source_to_kmer_count[seq_id] = number_of_kmers

# 3 pairwise reading
graph = rx.PyGraph()
nodes_indeces = graph.add_nodes_from(vertices)

def def_value(): return 0


batch_size = 1000000
batch_counter = 0
edges_tuples = []


print("[i] constructing graph")
with open(pairwise_file, 'r') as pairwise_tsv:
    next(pairwise_tsv)  # skip header

    for row in pairwise_tsv:
        row = row.strip().split()
        seq1 = int(row[1])
        seq2 = int(row[2])
        shared_kmers = int(row[3])
        min_seq = float(min(source_to_kmer_count[seq1], source_to_kmer_count[seq2]))
        containment = (shared_kmers / min_seq) * 100

        if containment < CONTAINMENT_THRESHOLD:
            continue
        
        if batch_counter < batch_size:
            batch_counter += 1
            edges_tuples.append((seq1 -1, seq2 - 1, containment))
        else:
            graph.add_edges_from(edges_tuples)
            batch_counter = 0
            edges_tuples.clear()

    else:
        if len(edges_tuples):
            graph.add_edges_from(edges_tuples)


no_connected_components = rx.number_connected_components(graph)
print(f"connected components: {no_connected_components}")


connected_components = rx.connected_components(graph)

with open(index_prefix + f"retworkx_{CONTAINMENT_THRESHOLD}.txt", 'w') as CLUSTERS:
    for component in connected_components:
        named_component = [id_to_name[node + 1] for node in component]
        CLUSTERS.write(','.join(named_component) + '\n')
        

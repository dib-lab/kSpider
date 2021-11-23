import sys
import os
import re

if len(sys.argv) < 3:
    sys.exit("run: python clusters_convert.py <clusters_tsv> <sraTable_csv>")

clusters_file = sys.argv[1]
sraRunTable = sys.argv[2]

# 0,16
runID_to_biome = dict()

# 0, 21
runID_to_geolocation = dict()

# 0, 3
runID_to_biosample = dict()

# combined
run_ID_combined = dict()


# with open(sraRunTable) as RUN_TABLE:
#     next(RUN_TABLE) # skip header
#     for line in RUN_TABLE:
#         line = line.strip()
#         line = re.split(r',(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)', line)
#         run_id = line[0]
#         biosample = line[3]
#         biome = line[16]
#         geolocation = line[21]
#         run_ID_combined[run_id] = f"{biosample}_{biome}_{geolocation}"
#         runID_to_geolocation[run_id] = geolocation
#         runID_to_biome[run_id] = biome


# output_biome = os.path.join(os.path.dirname(clusters_file), f"biome_{os.path.basename(clusters_file)}")
# output_geolocation = os.path.join(os.path.dirname(clusters_file), f"geolocation{os.path.basename(clusters_file)}")


output_combined = os.path.join(os.path.dirname(clusters_file), f"insights_{os.path.basename(clusters_file)}")
with open(clusters_file) as CLUSTERS, open(output_combined, 'w') as COMBINED:
    header = next(COMBINED)
    COMBINED.write(header)
    for line in CLUSTERS:
        line = line.strip().split('\t')
        cluster_id = line[0]
        all_clusters = line[1].split('|')
        new_clusters = [run_ID_combined[_cluster] for _cluster in all_clusters]
        new_clusters = '|'.join(new_clusters)
        COMBINED.write(f"{cluster_id}\t{new_clusters}\n")
        
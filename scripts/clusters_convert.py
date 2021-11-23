import sys
import os
import re

if len(sys.argv) < 3:
    sys.exit("run: python clusters_convert.py <clusters_tsv> <sraTable_csv>")

clusters_file = sys.argv[1]
sraRunTable = sys.argv[2]

# 0,16
runID_to_biome = dict()

# 9, 21
runID_to_geolocation = dict()

with open(sraRunTable) as RUN_TABLE:
    next(RUN_TABLE) # skip header
    for line in RUN_TABLE:
        line = line.strip()
        line = re.split(r',(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)', line)
        run_id = line[0]
        biome = line[16]
        geolocation = line[21]
        runID_to_geolocation[run_id] = geolocation
        runID_to_biome[run_id] = biome


output_biome = os.path.join(os.path.dirname(clusters_file), f"biome_{os.path.basename(clusters_file)}")
output_geolocation = os.path.join(os.path.dirname(clusters_file), f"geolocation{os.path.basename(clusters_file)}")



with open(clusters_file) as CLUSTERS, open(output_biome, 'w') as BIOME, open(output_geolocation, 'w') as GEO:
    header = next(CLUSTERS)
    BIOME.write(header)
    GEO.write(header)
    for line in CLUSTERS:
        line = line.strip().split('\t')
        cluster_id = line[0]
        all_clusters = line[1].split('|')
        biome_new_clusters = [runID_to_biome[_cluster] for _cluster in all_clusters]
        geo_new_clusters = [runID_to_geolocation[_cluster] for _cluster in all_clusters]
        biome_new_clusters = '|'.join(biome_new_clusters)
        geo_new_clusters = '|'.join(geo_new_clusters)
        BIOME.write(f"{cluster_id}\t{biome_new_clusters}\n")
        GEO.write(f"{cluster_id}\t{geo_new_clusters}\n")
        
        


    
    


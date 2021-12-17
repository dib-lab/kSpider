import sys
import os
import re

if len(sys.argv) < 4:
    sys.exit("run: python pairwise_ids_conversoin.py <pairwise_tsv> <namesMap_file> <sraTable_csv>")

pairwise_file = sys.argv[1]
namesMap = sys.argv[2]
sraRunTable = sys.argv[3]

# 0,16
runID_to_biome = dict()

# 0, 21
runID_to_geolocation = dict()

# 0, 3
runID_to_biosample = dict()

# OUTPUT
run_ID_OUTPUT = dict()


with open(sraRunTable) as RUN_TABLE:
    next(RUN_TABLE) # skip header
    for line in RUN_TABLE:
        line = line.strip()
        line = re.split(r',(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)', line)
        run_id = line[0]
        biosample = line[3]
        biome = line[16]
        geolocation = line[21]
        run_ID_OUTPUT[run_id] = f"{biosample}_{biome}_{geolocation}"
        runID_to_geolocation[run_id] = geolocation
        runID_to_biome[run_id] = biome


# output_biome = os.path.join(os.path.dirname(pairwise_file), f"biome_{os.path.basename(pairwise_file)}")
# output_geolocation = os.path.join(os.path.dirname(pairwise_file), f"geolocation{os.path.basename(pairwise_file)}")

"""
parse namesMap
"""

namesMap_dict = dict()
with open(namesMap) as NAMES:
    next(NAMES)
    for line in NAMES:
        line = line.strip().split()
        _id = line[0]
        _name = line[1]
        namesMap_dict[_id] = _name
    


converted_output = os.path.join(os.path.dirname(pairwise_file), f"converted_{os.path.basename(pairwise_file)}")

with open(pairwise_file) as CLUSTERS, open(converted_output, 'w') as OUTPUT:
    header = next(CLUSTERS) # skip headers line
    OUTPUT.write(header)
    for line in CLUSTERS:
        line = line.strip().split('\t')
        sample1 = namesMap_dict[line[1]] + '_' + run_ID_OUTPUT[namesMap_dict[line[1]]]
        sample2_id = namesMap_dict[line[2]] + '_' + run_ID_OUTPUT[namesMap_dict[line[2]]]
        shared_no_kmers = line[3]
        OUTPUT.write(f"{line[0]}\t{sample1}\t{sample2_id}\t{shared_no_kmers}\n")
        
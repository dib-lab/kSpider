import pickle
import sys
import os
from math import isclose

"""
LOAD GOLDEN FILES
"""

def count_lines_in_file(file_name):
    with open(file_name) as R:
        return sum(1 for _ in R)

with open('golden_pairwise.pickle', 'rb') as handle:
    golden_pairwise = pickle.load(handle)
    
with open('golden_sig_to_len.pickle', 'rb') as handle:
    golden_sig_to_len = pickle.load(handle)

with open("golden_min_containments.pickle", "rb") as handle:
    golden_min_containments = pickle.load(handle)

with open("golden_max_containments.pickle", "rb") as handle:
    golden_max_containments = pickle.load(handle)

with open("golden_avg_containments.pickle", "rb") as handle:
    golden_avg_containments = pickle.load(handle)

"""
LOAD kSpider files
"""

# Loading namesmap

sig_to_id = {}
id_to_sig = {}
with open('sigs.namesMap') as R:
    next(R)
    for line in R:
        line = line.strip().split()
        name = line[0]        
        sig_id = line[1]
        sig_to_id[sig_id] = name
        id_to_sig[name] = sig_id

# loading lengths
sig_to_len = {}
print("validating lengths")
with open('sigs_kSpider_seqToKmersNo.tsv') as R:
    next(R)
    for line in R:
        line = line.strip().split('\t')
        sig_id = id_to_sig[line[1]]
        sig_len = int(line[2])
        sig_to_len[sig_id] = sig_len
        
# loading containments
# source_1	source_2	shared_kmers	min_containment	avg_containment	max_containment
min_containments = {}
avg_containments = {}
max_contaiments = {}
pairwise = {}
kSpider_compare_sets = {}

if count_lines_in_file('sigs_kSpider_pairwise.tsv') < 100:
    print("ERROR! pairwise file is empty", file = sys.stderr)
    sys.exit(1)

with open('sigs_kSpider_pairwise.tsv') as R:
    next(R)
    for line in R:
        line = line.strip().split('\t')
        sig1 = id_to_sig[line[0]]
        sig2 = id_to_sig[line[1]]
        key = (sig1, sig2)
        shared_kmers = int(line[2])
        min_containment = float(line[3][:5])
        avg_containment = float(line[4][:5])
        max_containment = float(line[5][:5])
        min_containments[key] = min_containment
        avg_containments[key] = avg_containment
        max_contaiments[key] = max_containment
        pairwise[key] = shared_kmers

print("Loading kSpider's files completed")

"""
VALIDATING KMER COUNTS
"""

for sig, sig_len in sig_to_len.items():
    if sig_len != golden_sig_to_len[sig]:
        print(f"ERROR! {sig} lengths missmatch", file = sys.stderr)
else:
    print("kmer counts passed")

"""
VALIDATING PAIRWISE golden_pairwise
"""

shared_kmers_mismatches = 0
for sig_pair, shared_kmers in pairwise.items():
    if shared_kmers != golden_pairwise[sig_pair]:
        shared_kmers_mismatches += 1

if shared_kmers_mismatches:
    print(f"ERROR! {shared_kmers_mismatches} shared_kmers missmatch", file = sys.stderr)
else:
    print("pairwise validation passed")

"""
VALIDATING containments
"""

lengths_missmatch = 0

# validate lengths

for sig_id in sig_to_len:
    if sig_to_len[sig_id] != golden_sig_to_len[sig_id]:
        lengths_missmatch += 1

if lengths_missmatch:
    print(f"ERROR! {lengths_missmatch} lengths missmatch")

min_containments_mismatches = 0
avg_containments_mismatches = 0
max_containments_mismatches = 0
error_sigs = set()

# validate containments
for key in min_containments:
    if round(min_containments[key], 1) != round(golden_min_containments[key], 1):
        print(key, min_containments[key], golden_min_containments[key])
        error_sigs.add(key[0])
        error_sigs.add(key[1])
        min_containments_mismatches+=1    
    if round(avg_containments[key], 1) != round(golden_avg_containments[key], 1):
        error_sigs.add(key[0])
        error_sigs.add(key[1])
        avg_containments_mismatches+=1
    if round(max_contaiments[key], 1) != round(golden_max_containments[key], 1):
        max_containments_mismatches+=1
        error_sigs.add(key[0])
        error_sigs.add(key[1])
        
if min_containments_mismatches:
    print(f"ERROR! {min_containments_mismatches} min containments missmatch", file = sys.stderr)
if avg_containments_mismatches:
    print(f"ERROR! {avg_containments_mismatches} avg containments missmatch", file = sys.stderr)
if max_containments_mismatches:
    print(f"ERROR! {max_containments_mismatches} max containments missmatch", file = sys.stderr)

SUBSET_SIGS = False

if SUBSET_SIGS:
    # delete directory if exists
    if os.path.exists("subset_sigs"):
        os.rmdir("subset_sigs")
    
    os.mkdir("subset_sigs/")
        
    for sig in error_sigs:
        sig_file = f"sigs/{sig}.sig"
        # copy sig file to subset_sigs
        os.system(f"cp {sig_file} subset_sigs")

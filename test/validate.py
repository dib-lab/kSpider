import pickle
import sys

"""
LOAD GOLDEN FILES
"""

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

for sig_pair, shared_kmers in pairwise.items():
    if shared_kmers != golden_pairwise[sig_pair]:
        print(f"ERROR! {sig_pair} shared_kmers missmatch", file = sys.stderr)
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

# validate containments
for key in min_containments:
    if min_containments[key] != golden_min_containments[key]:
        min_containments_mismatches+=1    
    if avg_containments[key] != golden_avg_containments[key]:
        avg_containments_mismatches+=1
    if max_contaiments[key] != golden_max_containments[key]:
        max_containments_mismatches+=1
        
if min_containments_mismatches:
    print(f"ERROR! {min_containments_mismatches} min containments missmatch", file = sys.stderr)
if avg_containments_mismatches:
    print(f"ERROR! {avg_containments_mismatches} avg containments missmatch", file = sys.stderr)
if max_containments_mismatches:
    print(f"ERROR! {max_containments_mismatches} max containments missmatch", file = sys.stderr)

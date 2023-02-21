import json
from glob import glob
import os
from tqdm import tqdm
import pickle

def load_mins(json_file):
    with open(json_file) as f:
        data = json.load(f)
    return set(data[0]["signatures"][0]["mins"])

"""
LOADING
"""

# Load the signatures
sigs_list = glob("/home/mabuelanin/dib-dev/kSpider_dev/test/sigs/*sig")
sig_to_mins = {os.path.basename(file_name).replace('.sig',''): load_mins(file_name) for file_name in sigs_list}

"""
LENGTHS
"""
sig_to_len = {sig: len(sig_to_mins[sig]) for sig in sig_to_mins}

# Write lengths
with open('golden_sig_to_len.tsv', 'w') as f:
    for sig, sig_len in sig_to_len.items():
        f.write(f"{sig}\t{sig_len}\n")

with open('golden_sig_to_len.pickle', 'wb') as handle:
    pickle.dump(sig_to_len, handle, protocol=pickle.HIGHEST_PROTOCOL)


"""
SHARED_KMERS
"""

print("calculating shared kmers")
shared_kmers = {}
for sig1 in tqdm(sig_to_mins):
    for sig2 in sig_to_mins:
        if sig1 != sig2:
            key1 = (sig1, sig2)
            key2 = (sig1, sig2)
            
            common = len(sig_to_mins[sig1].intersection(sig_to_mins[sig2]))
            if common:
                shared_kmers[key1] = common
                shared_kmers[key2] = common


with open('golden_pairwise.pickle', 'wb') as handle:
    pickle.dump(shared_kmers, handle, protocol=pickle.HIGHEST_PROTOCOL)

"""
CONTAINMENT
"""

print("calculating containments")

min_containments = {}
max_containments = {}
avg_containments = {}


for sig1 in tqdm(sig_to_mins):
    for sig2 in sig_to_mins:
        if sig1 != sig2:
            key1 = (sig1, sig2)
            key2 = (sig1, sig2)
            
            common = shared_kmers[key1]
            if common:                
                len_sig1 = sig_to_len[sig1]
                len_sig2 = sig_to_len[sig2]
                max_containment = float(common) / min(len_sig1, len_sig2)
                min_containment = float(common) / max(len_sig1, len_sig2)
                avg_containment = (max_containment + min_containment) / 2
                
                min_containment = float(f"{min_containment:.3f}")
                avg_containment = float(f"{avg_containment:.3f}")
                max_containment = float(f"{max_containment:.3f}")
                
                min_containments[key1] = min_containment
                min_containments[key2] = min_containment
                avg_containments[key1] = avg_containment
                avg_containments[key2] = avg_containment
                max_containments[key1] = max_containment
                max_containments[key2] = max_containment

with open('golden_min_containments.pickle', 'wb') as handle:
    pickle.dump(min_containments, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open('golden_max_containments.pickle', 'wb') as handle:
    pickle.dump(max_containments, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open('golden_avg_containments.pickle', 'wb') as handle:
    pickle.dump(avg_containments, handle, protocol=pickle.HIGHEST_PROTOCOL)
    

# max_containments to tsv file
with open('golden_max_containments.tsv', 'w') as f:
    for key, value in max_containments.items():
        f.write(f"{key[0]}\t{key[1]}\t{value}")
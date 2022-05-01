import json
import glob
import gzip
import io
import codecs
import sys
import os


def get_hashes(filename, kSize):

    def is_gz_file(filepath):
        with open(filepath, 'rb') as test_f:
            return test_f.read(2) == b'\x1f\x8b'

    def select(sigs):
        for sig in sigs:
            if sig["ksize"] == kSize:
                return set(sig["mins"])

    if is_gz_file(filename):
        with gzip.open(filename, 'r') as fin:
            signatures = json.loads(
                fin.read().decode('utf-8'))[0]["signatures"]
            return select(signatures)

    else:
        with open(filename, 'r') as fin:
            signatures = json.loads(
                fin.read().decode('utf-8'))[0]["signatures"]
            return select(signatures)


if len(sys.argv) < 3:
    sys.exit("usage: Python <dir:sigs_dir> <int:kSize>")


file_reader = None
json_str = ""
file_to_hashes = {}
sigs_dir = sys.argv[1]
kSize = int(sys.argv[2])
namesmap = {}
containment = {}


for i, sig_file in enumerate(glob.glob(sigs_dir + "/*")):
    namesmap[i] = os.path.basename(sig_file)
    file_to_hashes[i] = get_hashes(sig_file, kSize)

all_ids = list(namesmap.keys())
pairs = [(a, b) for idx, a in enumerate(all_ids) for b in all_ids[idx + 1:]]


for pair in pairs:
    shared_kmers = len(
        file_to_hashes[pair[0]].intersection(file_to_hashes[pair[1]]))
    if shared_kmers:
        containment[pair] = shared_kmers

print("seq1\tseq2\tshared_kmers")
for (grp_id_1, grp_id_2), shared_kmers in containment.items():
    print(f"{grp_id_1}[{namesmap[grp_id_1]}]\t{grp_id_2}[{namesmap[grp_id_2]}]\t{shared_kmers}")

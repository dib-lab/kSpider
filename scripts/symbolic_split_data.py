"""[summary]
You have so many sequence files and want to split in subdirs for parallelization?, Then, this script for you.

[ ] Support paired-end
[ ] Support single end.
"""

from glob import glob
import os
import sys


if len(sys.argv) < 3 or sys.argv[1] != "-c":
    sys.exit("usage: python symbolic_split_data.py -c <chunks> <seqs_dir>")


chunks = int(sys.argv[2])
data = sorted(list(glob(os.path.join(sys.argv[3], '*'))))


BOOL_PAIRED_END = False

## verify it's paired_end
_f1 = os.path.basename(data[0])
_f2 = os.path.basename(data[1])

_d = []
for i in range(len(_f1)): 
    if _f1[i] != _f2[i]:
        _d.append(_f1[i])
        _d.append(_f2[i])
        
if "1" in _d and "2" in _d:
    BOOL_PAIRED_END = True
else:
    print("Not Implemented!")

## Start processing

if BOOL_PAIRED_END:
    read_pairs = list()
    read_groups = list()
    for i in range(0, len(data), 2):
        read_pairs.append((data[i], data[i+1]))
    
    n = chunks
    read_groups = [read_pairs[i * chunks:(i + 1) * chunks] for i in range((len(read_pairs) + chunks - 1) // chunks )] 
    subdir_prefix = os.path.dirname(sys.argv[3])
    
    for i in range(len(read_groups)):
        pwd = os.path.abspath(os.getcwd())
        dir_name = os.path.join(pwd, f"{subdir_prefix}_{i+1}")
        os.mkdir(dir_name)
        for r1, r2 in read_groups[i]:
            r1 = os.path.abspath(r1)
            r2 = os.path.abspath(r2)
            ln_r1 = os.path.join(dir_name, os.path.basename(r1))
            ln_r2 = os.path.join(dir_name, os.path.basename(r2))
                   
            if os.path.islink(ln_r1):
                os.unlink(ln_r1)
            os.symlink(r1, ln_r1)
            os.symlink(r2, ln_r2)

    print("splitting done successfully")     

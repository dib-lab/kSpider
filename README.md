# kSpider2 | C++ verison of the kSpider software tool.

## Installation

### Clone

```bash
git clone https://github.com/mr-eyes/kCluster2.git
cd kCluster2
git submodule update --init --recursive
```

### Install

```bash
mkdir build && cd build/
cmake .. && make
cd ..
```

## Usage

```bash

MIN_Q=5
MAX_Q=25
STEP_Q=5

INDEX_PREFIX=idx_seq

./kCluster_pairwise --min-q=${MIN_Q} --max-q=${MAX_Q} --step-q=${STEP_Q} --idx=${INDEX_PREFIX}

```

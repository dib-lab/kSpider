# kSpider2 | C++ verison of the kSpider software tool.

## Installation

### Clone

```bash
git clone https://github.com/mr-eyes/kSpider2.git
cd kSpider2
git submodule update --init --recursive
```

### Install

```bash
PROJECT=$(pwd)
cd lib/kProcessor
mkdir build && cd build && cmake .. && make
cd $PROJECT
mkdir build && cd build/
cmake .. && make
cd ..
```

## Usage

```bash

INDEX_PREFIX=idx_seq

./kSpider_pairwise ${INDEX_PREFIX}

```

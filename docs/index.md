# kSpider

## Description

First, it creates an index using kProcessor for the source sequences. Second, it constructs a pairwise containment matrix through a single iteration over the index. Finally, it builds a graph from the pairwise matrix and applies a connected-components graph algorithm to extract the clusters with a user-defined containment threshold.

## Understanding the Indexing

The first step in clustering sequences is to index them. Indexing can be done in multiple flavors, depending on your use case. For instance, you have a transcriptome, and you want to perform clustering on both transcript and gene levels. Both indexing processes share the same FASTA file. To discriminate between the transcript and gene levels, we will need to create another file called "names file". This **tab-separated** file contains two columns. The first column must contain the same FASTA header without the '>', and the second column contains the user-defined group name corresponding to the FASTA header. So, in the transcript-level indexing, the second column will be the transcript ID. However, in the gene-level indexing, the second column will contain the gene ID.

**FASTA FILE**
```txt
>TR1|GENE1
AGCTAGACTACTACTACGACTAGGAGACTCA
>TR2|GENE1
AGCTAGCTATACGATCGTACATCGAATCAAC
>TR3|GENE2
GTGATGACTGATCGATAATGCTAGCCACATG
>TR4|GENE2
ATCACGTAGCTATTATCGATCGACTACACAC
```

**Transcript-level names file**
```txt
TR1|GENE1	TR1
TR2|GENE1	TR2
TR3|GENE2	TR3
TR4|GENE2	TR4
```

**Gene-level names file**
```txt
TR1|GENE1	GENE1
TR2|GENE1	GENE1
TR3|GENE2	GENE2
TR4|GENE2	GENE2
```

---

## Installation

### Easy installation

```bash
pip install kSpider
```

### Manual Installation / Development

```bash

Clone and build
git clone https://github.com/mr-eyes/kSpider.git
cd kSpider
cmake -Bbuild
cmake --build build
bash build_wrapper.sh
```

---

## Indexing

kSpider supports 3 main indexing mode.
1. Kmers mode, where you can select the kmer size.
2. Skipmers mode, where you need to provide m,n,k values. Read more at https://www.biorxiv.org/content/10.1101/179960v2
3. Protein mode, where you can set a kSize up to 11 amino acid in the default mode, or 19 amino acid in the dayhoff mode.



## Usage Example

**Download human protein-coding sequences**
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.pc_transcripts.fa.gz
```

**Generate names files**
```bash
# Generate for genes indexing
zcat gencode.v38.pc_transcripts.fa | grep ">" | cut -c2- |  awk -F'|' '{print $0"\t"$2}' > group_by_genes.names
```

**Index kmers by genes**
```bash
kSpider index_kmers -f gencode.v38.pc_transcripts.fa.gz -n group_by_genes.names -k 25 -o idx_genes
```

**Generate pairwise containment matrix for the two indeces**
```bash
kSpider pairwise -i idx_genes
```

**Perform clustering at %90 containment threshold**
```bash
kSpider cluster -i idx_genes -c 0.95
```
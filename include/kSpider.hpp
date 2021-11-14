#include <iostream>
#include <cstdint>
#include "argh.h"
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"


namespace kSpider{

void pairwise(string index_prefix);
void index_kmers(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
void index_kmers_nonCanonical(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
void index_skipmers(int m, int n, int k, string fasta_file, string names_file, int chunk_size, string index_prefix);
void index_protein(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
void index_dayhoff(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);

};
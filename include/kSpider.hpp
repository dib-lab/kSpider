#include <iostream>
#include <cstdint>
#include "argh.h"
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include <queue>

namespace kSpider {

    void pairwise(string index_prefix);
    void index_kmers(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    void index_kmers_nonCanonical(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    void index_skipmers(int m, int n, int k, string fasta_file, string names_file, int chunk_size, string index_prefix);
    void index_protein(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    void index_dayhoff(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);

    void index_datasets(string kfs_dir);

    void paired_end_to_kDataFrame(string r1_file_name, string r2_file_name, int kSize, int chunk_size);
    void single_end_to_kDataFrame(string r1_file_name, int kSize, int chunk_size);
    void protein_to_kDataFrame(string r1_file_name, int kSize, int chunk_size, bool is_dayhoff, string output_prefix);

};
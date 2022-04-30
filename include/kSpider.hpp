#include <iostream>
#include <cstdint>
#include "argh.h"
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include <queue>
#include <string>

using namespace std;

namespace kSpider {

    void pairwise(string index_prefix);
    void index_kmers(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    void index_kmers_nonCanonical(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    void index_skipmers(int m, int n, int k, string fasta_file, string names_file, int chunk_size, string index_prefix);
    void index_protein(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    void index_dayhoff(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    void index_datasets(string kfs_dir);

    // The new one (index priotiyQueue)
    void datasets_indexing(string kfs_dir);

    void sourmash_index_kp1(string sigs_dir, int selective_kSize);
    void sourmash_index_kp1_fast(string sigs_dir, int selective_kSize);
    

    void paired_end_to_kDataFrame(string r1_file_name, string r2_file_name, int kSize, int chunk_size, int downsampling_ratio, bool remove_singletones, bool ondisk);
    void single_end_to_kDataFrame(string r1_file_name, int kSize, int chunk_size, int downsampling_ratio, bool remove_singletones, bool ondisk);
    void protein_to_kDataFrame(string r1_file_name, int kSize, int chunk_size, bool is_dayhoff, string output_prefix, int downsampling_ratio, bool ondisk);

}
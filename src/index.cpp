#include "kSpider.hpp"
#include <iostream>
#include <cstdint>
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include "algorithms.hpp"

namespace kSpider {

    void index_kmers(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix) {
        auto* KF = new kDataFrameMQF(KMERS, integer_hasher, { {"kSize", kSize} });
        auto* ckf = kProcessor::index(KF, fasta_file, chunk_size, names_file);
        ckf->save(index_prefix);
    }

    void index_kmers_nonCanonical(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix) {
        auto* KF = new kDataFrameMQF(KMERS, nonCanonicalInteger_Hasher, { {"kSize", kSize} });
        auto* ckf = kProcessor::index(KF, fasta_file, chunk_size, names_file);
        ckf->save(index_prefix);
    }

    void index_skipmers(int m, int n, int k, string fasta_file, string names_file, int chunk_size, string index_prefix) {
        auto* KF = new kDataFramePHMAP(SKIPMERS, integer_hasher, { {"m", m}, {"n", n}, {"k", k} });
        auto* ckf = kProcessor::index(KF, fasta_file, chunk_size, names_file);
        ckf->save(index_prefix);
    }

    void index_protein(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix) {
        auto* KF = new kDataFramePHMAP(PROTEIN, protein_hasher, { {"kSize", kSize} });
        auto* ckf = kProcessor::index(KF, fasta_file, chunk_size, names_file);
        ckf->save(index_prefix);
    }

    void index_dayhoff(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix) {
        auto* KF = new kDataFramePHMAP(PROTEIN, proteinDayhoff_hasher, { {"kSize", kSize} });
        auto* ckf = kProcessor::index(KF, fasta_file, chunk_size, names_file);
        ckf->save(index_prefix);
    }

}
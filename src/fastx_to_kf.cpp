
#include <iostream>
#include <kDataFrame.hpp>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithms.hpp>
#include "assert.h"
#include <cstring>
#include "kSpider.hpp"
#include <cstdint>
#include <math.h>

namespace kSpider {

    void paired_end_to_kDataFrame(string r1_file_name, string r2_file_name, int kSize, int chunk_size, int downsampling_ratio, bool remove_singletones) {

        string PE_1_reads_file = r1_file_name;
        string PE_2_reads_file = r2_file_name;

        std::string base_filename = PE_1_reads_file.substr(PE_1_reads_file.find_last_of("/\\") + 1);
        base_filename = base_filename.substr(0, base_filename.find('_'));

        kmerDecoder* READ_1_KMERS = kmerDecoder::getInstance(r1_file_name, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });
        kmerDecoder* READ_2_KMERS = kmerDecoder::getInstance(r2_file_name, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });

        auto* kf = new kDataFramePHMAP(KMERS, mumur_hasher, { {"kSize", kSize} });

        int Reads_chunks_counter = 0;
        uint64_t max_hash = UINT64_MAX / (uint64_t)downsampling_ratio;
        uint64_t total_kmers = 0;
        uint64_t inserted_kmers = 0;


        while (!READ_1_KMERS->end() && !READ_2_KMERS->end()) {


            READ_1_KMERS->next_chunk();
            READ_2_KMERS->next_chunk();


            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1 = READ_1_KMERS->getKmers()->begin();
            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq2 = READ_2_KMERS->getKmers()->begin();
            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1_end = READ_1_KMERS->getKmers()->end();
            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq2_end = READ_2_KMERS->getKmers()->end();

            while (seq1 != seq1_end && seq2 != seq2_end) {

                for (auto const kRow : seq1->second) {
                    if (kRow.hash < max_hash) {
                        kf->insert(kRow.hash);
                        total_kmers++;
                        inserted_kmers++;
                    }
                    else {
                        total_kmers++;
                        continue;
                    }
                }


                for (auto const kRow : seq2->second) {
                    if (kRow.hash < max_hash) {
                        kf->insert(kRow.hash);
                        total_kmers++;
                        inserted_kmers++;
                    }
                    else {
                        total_kmers++;
                        continue;
                    }
                }

                seq1++;
                seq2++;

            }

        }
        int removed_singletones = 0;
        if (remove_singletones) {
            auto* new_kf = new kDataFramePHMAP(kSize);
            auto it = kf->begin();
            while (it != kf->end()) {
                if (it.getCount() > 1) {
                    new_kf->insert(it.getHashedKmer(), it.getCount());
                    it++;
                }
                removed_singletones++;
                it++;
            }
            cout << "removed " << removed_singletones << " singletones." << endl;
            new_kf->save(base_filename);
            cout << "filename(" << base_filename << "): total(" << total_kmers << ") inserted(" << (inserted_kmers - removed_singletones) << ") << inserted_unique(" << new_kf->size() << ")" << endl;
        }
        else {
            kf->save(base_filename);
            cout << "filename(" << base_filename << "): total(" << total_kmers << ") inserted(" << inserted_kmers << ") << inserted_unique(" << kf->size() << ")" << endl;
        }
    }

    void single_end_to_kDataFrame(string r1_file_name, int kSize, int chunk_size, int downsampling_ratio, bool remove_singletones) {

        string PE_1_reads_file = r1_file_name;

        std::string base_filename = PE_1_reads_file.substr(PE_1_reads_file.find_last_of("/\\") + 1);

        kmerDecoder* READ_1_KMERS = kmerDecoder::getInstance(r1_file_name, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });

        auto* kf = new kDataFramePHMAP(KMERS, mumur_hasher, { {"kSize", kSize} });

        int Reads_chunks_counter = 0;
        uint64_t max_hash = UINT64_MAX / (uint64_t)downsampling_ratio;
        uint64_t total_kmers = 0;
        uint64_t inserted_kmers = 0;


        while (!READ_1_KMERS->end()) {


            READ_1_KMERS->next_chunk();

            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1 = READ_1_KMERS->getKmers()->begin();
            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1_end = READ_1_KMERS->getKmers()->end();

            while (seq1 != seq1_end) {

                for (auto const kRow : seq1->second) {
                    if (kRow.hash < max_hash) {
                        kf->insert(kRow.hash);
                        total_kmers++;
                        inserted_kmers++;
                    }
                    else {
                        total_kmers++;
                        continue;
                    }
                }

                seq1++;
            }

        }
        int removed_singletones = 0;
        if (remove_singletones) {
            auto* new_kf = new kDataFramePHMAP(kSize);
            auto it = kf->begin();
            while (it != kf->end()) {
                if (it.getCount() > 1) {
                    new_kf->insert(it.getHashedKmer(), it.getCount());
                    it++;
                }
                removed_singletones++;
                it++;
            }
            cout << "removed " << removed_singletones << " singletones." << endl;
            new_kf->save(base_filename);
            cout << "filename(" << base_filename << "): total(" << total_kmers << ") inserted(" << (inserted_kmers - removed_singletones) << ") << inserted_unique(" << new_kf->size() << ")" << endl;
        }
        else {
            kf->save(base_filename);
            cout << "filename(" << base_filename << "): total(" << total_kmers << ") inserted(" << inserted_kmers << ") << inserted_unique(" << kf->size() << ")" << endl;
        }
    }


    void protein_to_kDataFrame(string r1_file_name, int kSize, int chunk_size, bool is_dayhoff, string output_prefix, int downsampling_ratio) {

        string PE_1_reads_file = r1_file_name;

        hashingModes hasher_type = protein_hasher;
        if (is_dayhoff) hasher_type = proteinDayhoff_hasher;

        kmerDecoder* KD = kmerDecoder::getInstance(r1_file_name, chunk_size, PROTEIN, hasher_type, { {"kSize", kSize} });
        kDataFramePHMAP* kf = new kDataFramePHMAP(PROTEIN, hasher_type, { {"kSize", kSize} });

        int hashing_kmer_kSize = (int)((kSize * 5) / 2);


        auto* INT_HASHER = new IntegerHasher(hashing_kmer_kSize);
        uint64_t max_real_hash = INT_HASHER->hash(pow(2, hashing_kmer_kSize));

        uint64_t max_hash = max_real_hash / downsampling_ratio;

        if (downsampling_ratio == 1) max_hash = UINT64_MAX;

        uint64_t total_kmers = 0;
        uint64_t inserted_kmers = 0;

        while (!KD->end()) {
            KD->next_chunk();

            for (const auto& seq : *KD->getKmers()) {
                for (const auto& kmer : seq.second) {
                    // downsampling
                    uint64_t kmer_hash = INT_HASHER->hash(kmer.hash);
                    // cout << kmer_hash << endl;
                    if (kmer_hash < max_hash) {
                        kf->insert(kmer.hash); // Insert the 5-bit representation not the hash val.
                        total_kmers++;
                        inserted_kmers++;
                    }
                    else {
                        total_kmers++;
                        continue;
                    }
                }
            }
        }

        cout << "filename(" << output_prefix << "): total(" << total_kmers << ") inserted(" << inserted_kmers << ")" << endl;
        kf->save(output_prefix);
    }

}

#include <iostream>
#include <kDataFrame.hpp>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithms.hpp>
#include "assert.h"
#include <cstring>
#include "kSpider.hpp"

namespace kSpider {

    void paired_end_to_kDataFrame(string r1_file_name, string r2_file_name, int kSize, int chunk_size) {

        string PE_1_reads_file = r1_file_name;
        string PE_2_reads_file = r2_file_name;

        std::string base_filename = PE_1_reads_file.substr(PE_1_reads_file.find_last_of("/\\") + 1);
        base_filename = base_filename.substr(0, base_filename.find('_'));

        kmerDecoder* READ_1_KMERS = kmerDecoder::getInstance(r1_file_name, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });
        kmerDecoder* READ_2_KMERS = kmerDecoder::getInstance(r2_file_name, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });

        int Reads_chunks_counter = 0;


        while (!READ_1_KMERS->end() && !READ_2_KMERS->end()) {

            cerr << "processing chunk: (" << ++Reads_chunks_counter << "...";
            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

            READ_1_KMERS->next_chunk();
            READ_2_KMERS->next_chunk();


            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1 = READ_1_KMERS->getKmers()->begin();
            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq2 = READ_2_KMERS->getKmers()->begin();
            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1_end = READ_1_KMERS->getKmers()->end();
            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq2_end = READ_2_KMERS->getKmers()->end();

            while (seq1 != seq1_end && seq2 != seq2_end) {

                auto* kf = new kDataFrameMQF(KMERS, mumur_hasher, { {"kSize", kSize} });

                for (auto const kRow : seq1->second) kf->insert(kRow.hash);
                for (auto const kRow : seq2->second) kf->insert(kRow.hash);

                kf->save(base_filename);

                seq1++;
                seq2++;

                delete kf;
            }


            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
            long hr = milli / 3600000;
            milli = milli - 3600000 * hr;
            long min = milli / 60000;
            milli = milli - 60000 * min;
            long sec = milli / 1000;
            milli = milli - 1000 * sec;
            cerr << "Done in: ";
            cerr << min << ":" << sec << ":" << milli << endl;

        }

    }

    void single_end_to_kDataFrame(string r1_file_name, int kSize, int chunk_size) {

        string PE_1_reads_file = r1_file_name;

        std::string base_filename = PE_1_reads_file.substr(PE_1_reads_file.find_last_of("/\\") + 1);
        base_filename = base_filename.substr(0, base_filename.find('_'));

        kmerDecoder* READ_1_KMERS = kmerDecoder::getInstance(r1_file_name, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });

        int Reads_chunks_counter = 0;


        while (!READ_1_KMERS->end()) {

            cerr << "processing chunk: (" << ++Reads_chunks_counter << "...";
            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

            READ_1_KMERS->next_chunk();


            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1 = READ_1_KMERS->getKmers()->begin();
            flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1_end = READ_1_KMERS->getKmers()->end();

            while (seq1 != seq1_end) {
                auto* kf = new kDataFrameMQF(KMERS, integer_hasher, { {"kSize", kSize} });
                for (auto const kRow : seq1->second) {
                    // downsampling
                    if (kRow.hash < 147573952589676412) kf->insert(kRow.hash);
                }
                kf->save(base_filename);
                seq1++;
            }


            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
            long hr = milli / 3600000;
            milli = milli - 3600000 * hr;
            long min = milli / 60000;
            milli = milli - 60000 * min;
            long sec = milli / 1000;
            milli = milli - 1000 * sec;
            cerr << "Done in: ";
            cerr << min << ":" << sec << ":" << milli << endl;

        }

    }


    void protein_to_kDataFrame(string r1_file_name, int kSize, int chunk_size, bool is_dayhoff, string output_prefix) {

        string PE_1_reads_file = r1_file_name;

        hashingModes hasher_type = protein_hasher;
        if (is_dayhoff) hasher_type = proteinDayhoff_hasher;

        kmerDecoder* KD = kmerDecoder::getInstance(r1_file_name, chunk_size, PROTEIN, hasher_type, { {"kSize", kSize} });
        kDataFramePHMAP * kf = new kDataFramePHMAP(PROTEIN, hasher_type, { {"kSize", kSize} });


        while (!KD->end()) {
            KD->next_chunk();

            for (const auto& seq : *KD->getKmers()) {
                for (const auto& kmer : seq.second) {
                    // Downsampling
                    kf->insert(kmer.hash);
                    // if (kmer.hash < 147573952589676412) kf->insert(kmer.hash); else continue;
                }
            }
        }

        kf->save(output_prefix);

    }

}
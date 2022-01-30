#include "kDataFrame.hpp"
#include "algorithms.hpp"
#include <cmath>
#include <bits/stdc++.h>

#define LOWEST_PERCENTILE 5
#define DOWNSAMPLING_RATIO 250

int main(int argc, char** argv) {


    string PE_1_reads_file = argv[1];
    string PE_2_reads_file = argv[2];
    int chunk_size = 1000;
    int kSize = 25;
    uint64_t DESIRED_NUM_KMERS = 100000000; // 10 MILLIONS

    uint64_t stats_total_kmers = 0;
    uint64_t stats_total_kmers_unique = 0;
    uint64_t stats_singletones = 0;
    uint64_t stats_after_percentile_unique = 0;
    uint64_t stats_final_inserted = 0;
    uint64_t stats_final_inserted_unique = 0;
    uint64_t stats_removed_percentile = 0;



    std::string base_filename = PE_1_reads_file.substr(PE_1_reads_file.find_last_of("/\\") + 1);
    base_filename = base_filename.substr(0, base_filename.find('_'));

    cout << "SAMPLE: " << base_filename << endl;


    kmerDecoder* READ_1_KMERS = kmerDecoder::getInstance(PE_1_reads_file, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });
    kmerDecoder* READ_2_KMERS = kmerDecoder::getInstance(PE_2_reads_file, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });

    flat_hash_map<uint64_t, uint64_t> kf;

    int Reads_chunks_counter = 0;
    uint64_t total_kmers = 0;
    uint64_t inserted_kmers = 0;

    // First iteration: Insert and count the kmers in the kDataFrame
    cout << "First iteration: Insert and count the kmers in the kDataFrame" << endl;
    while (!READ_1_KMERS->end() && !READ_2_KMERS->end()) {
        READ_1_KMERS->next_chunk();
        READ_2_KMERS->next_chunk();
        flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1 = READ_1_KMERS->getKmers()->begin();
        flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq2 = READ_2_KMERS->getKmers()->begin();
        flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1_end = READ_1_KMERS->getKmers()->end();
        flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq2_end = READ_2_KMERS->getKmers()->end();
        while (seq1 != seq1_end && seq2 != seq2_end) {

            for (auto const kRow : seq1->second) {
                kf[kRow.hash]++;
                stats_total_kmers++;
            }
            for (auto const kRow : seq2->second) {
                kf[kRow.hash]++;
                stats_total_kmers++;
            }

            seq1++;
            seq2++;
        }
    }

    stats_total_kmers_unique = kf.size();

    cout << "Second iteration: remove lowest %" << LOWEST_PERCENTILE << " kmer count percentile" << endl;

    flat_hash_map<uint64_t, uint64_t>::iterator kf_it = kf.begin();
    auto* kmerCounts = new(vector<uint64_t>);
    while (kf_it != kf.end()) kmerCounts->push_back(kf_it->second);
    sort(kmerCounts->begin(), kmerCounts->end());
    uint64_t _idx = (uint64_t)ceil((kf.size() * LOWEST_PERCENTILE / 100));
    uint64_t count_percentile_cutoff = kmerCounts->at(_idx);
    cout << "    calculated percentile cutoff kmercount=" << count_percentile_cutoff << endl;
    delete kmerCounts;


    flat_hash_map<uint64_t, uint64_t>::iterator kf_it = kf.begin();
    while (kf_it != kf.end()) {
        if (kf_it->second <= count_percentile_cutoff) {
            kf_it = kf.erase(kf_it);
            stats_removed_percentile++;
        }
        else {
            kf_it++;
        }
    }

    stats_after_percentile_unique = kf.size();

    cout << "Third iteration: downsampling to 1/" << DOWNSAMPLING_RATIO << endl;

    uint64_t max_real_hash = READ_1_KMERS->hasher->hash(pow(2, kSize));
    uint64_t max_hash = UINT64_MAX / (uint64_t)DOWNSAMPLING_RATIO;

    // Check if the sample size is less than the desired num of kmers
    if (kf.size() == 0) {

        // Move into kDataFrame
        auto* kf_final = new kDataFramePHMAP(KMERS, mumur_hasher, { {"kSize", kSize} });

        flat_hash_map<uint64_t, uint64_t>::iterator _kf_it = kf.begin();
        while (_kf_it != kf.end()) {
            kf_final->insert(_kf_it->first, _kf_it->second);
            _kf_it = kf.erase(_kf_it);
        }

        cout << "saving sample (" << base_filename << ")" << endl;
        kf_final->save("normalized_" + base_filename);
        cout << "sample(" << base_filename << ") has (" << kf_final->size() << ") kmers which is less than the desired (" << DESIRED_NUM_KMERS << ") kmers.";
        cout << "stats_total_kmers= " << stats_total_kmers << endl;
        cout << "stats_total_kmers_unique= " << stats_total_kmers_unique << endl;
        cout << "stats_singletones= " << stats_singletones << endl;
        cout << "stats_after_percentile_unique= " << stats_after_percentile_unique << endl;
        cout << "__________________________________" << endl;
        exit(0);
    }

    // Third Iteration: Calculaute the lowest %5 percentile


    // Third Iteration: get only the number of kmers we want

    int step = static_cast<double>(stats_after_percentile_unique) / DESIRED_NUM_KMERS;

    cout << "Third iteration: get only the number of kmers we want" << endl;
    auto* final_kf = new kDataFramePHMAP(KMERS, mumur_hasher, { {"kSize", kSize} });

    flat_hash_map<uint64_t, uint64_t>::iterator _kf_it = kf.begin();
    while (_kf_it != kf.end() && final_kf->size() < DESIRED_NUM_KMERS) {
        final_kf->insert(_kf_it->first, _kf_it->second);
        _kf_it = kf.erase(_kf_it);
        for (int i = 0; i < step; i++) _kf_it++;
    }

    stats_final_inserted_unique = final_kf->size();

    cout << "stats_total_kmers= " << stats_total_kmers << endl;
    cout << "stats_total_kmers_unique= " << stats_total_kmers_unique << endl;
    cout << "stats_singletones= " << stats_singletones << endl;
    cout << "stats_after_percentile_unique= " << stats_after_percentile_unique << endl;
    cout << "stats_final_inserted_unique= " << stats_final_inserted_unique << endl;

    final_kf->save("normalized_" + base_filename);
    cout << "__________________________________" << endl;


}
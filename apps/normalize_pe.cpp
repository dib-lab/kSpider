#include "kDataFrame.hpp"
#include "algorithms.hpp"



int main(int argc, char** argv) {


    string PE_1_reads_file = ""; // r1_file_name;
    string PE_2_reads_file = ""; // r2_file_name;
    int chunk_size = 1000;
    int kSize = 25;
    uint64_t DESIRED_NUM_KMERS = 100;

    uint64_t stats_total_kmers = 0;
    uint64_t stats_total_kmers_unique = 0;
    uint64_t stats_singletones = 0;
    uint64_t stats_after_singletone_unique = 0;
    uint64_t stats_final_inserted = 0;
    uint64_t stats_final_inserted_unique = 0;


    std::string base_filename = PE_1_reads_file.substr(PE_1_reads_file.find_last_of("/\\") + 1);
    base_filename = base_filename.substr(0, base_filename.find('_'));

    kmerDecoder* READ_1_KMERS = kmerDecoder::getInstance(PE_1_reads_file, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });
    kmerDecoder* READ_2_KMERS = kmerDecoder::getInstance(PE_2_reads_file, chunk_size, KMERS, mumur_hasher, { {"kSize", kSize} });

    auto* kf = new kDataFramePHMAP(KMERS, mumur_hasher, { {"kSize", kSize} });

    int Reads_chunks_counter = 0;
    uint64_t total_kmers = 0;
    uint64_t inserted_kmers = 0;

    // First iteration: Insert and count the kmers in the kDataFrame
    while (!READ_1_KMERS->end() && !READ_2_KMERS->end()) {
        READ_1_KMERS->next_chunk();
        READ_2_KMERS->next_chunk();
        flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1 = READ_1_KMERS->getKmers()->begin();
        flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq2 = READ_2_KMERS->getKmers()->begin();
        flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq1_end = READ_1_KMERS->getKmers()->end();
        flat_hash_map<std::string, std::vector<kmer_row>>::iterator seq2_end = READ_2_KMERS->getKmers()->end();
        while (seq1 != seq1_end && seq2 != seq2_end) {

            for (auto const kRow : seq1->second) {
                kf->insert(kRow.hash);
                stats_total_kmers++;
            }
            for (auto const kRow : seq2->second) {
                kf->insert(kRow.hash);
                stats_total_kmers++;
            }

            seq1++;
            seq2++;
        }
    }

    stats_total_kmers_unique = kf->size();


    // Second iteration: remove singletones
    auto* kf_singletones_free = new kDataFramePHMAP(KMERS, mumur_hasher, { {"kSize", kSize} });

    auto kf_it = kf->begin();
    while (kf_it != kf->end()) {
        if (kf_it.getCount() > 1) kf_singletones_free->insert(kf_it.getHashedKmer());
        else stats_singletones++;
        kf_it++;
    }

    stats_after_singletone_unique = kf_singletones_free->size();

    // Now delete the old kDataFrame
    delete kf;

    // Check if the sample size is less than the desired num of kmers
    if (kf_singletones_free->size() < DESIRED_NUM_KMERS) {
        cout << "saving sample (" << base_filename << ")" << endl;
        kf_singletones_free->save(base_filename);
        cout << "sample() has (" << kf_singletones_free->size() << ") kmers which is less than the desired (" << DESIRED_NUM_KMERS << ") kmers.";
        exit(0);
    }

    // Third Iteration: get only the number of kmers we want

    int step = static_cast<double>(stats_after_singletone_unique) / DESIRED_NUM_KMERS;

    cout << "Third iteration: get only the number of kmers we want" << endl;
    auto* final_kf = new kDataFramePHMAP(KMERS, mumur_hasher, { {"kSize", kSize} });
    auto kf_singletones_free_it = kf->begin();

    while (kf_singletones_free_it != kf->end() && final_kf->size() < DESIRED_NUM_KMERS) {
        final_kf->insert(kf_singletones_free_it.getHashedKmer());
        for (int i = 0; i < step; i++) kf_singletones_free_it++;
    }

    stats_final_inserted_unique = final_kf->size();

    cout << "stats_total_kmers= " << stats_total_kmers << endl;
    cout << "stats_total_kmers_unique= " << stats_total_kmers_unique << endl;
    cout << "stats_singletones= " << stats_singletones << endl;
    cout << "stats_after_singletone_unique= " << stats_after_singletone_unique << endl;
    cout << "stats_final_inserted_unique= " << stats_final_inserted_unique << endl;

}
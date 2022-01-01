#include "kDataFrame.hpp"
#include "algorithms.hpp"
#include <iostream>
#include <fstream>
using namespace std;

#define LOWEST_PERCENTILE 5

int main(int args, char ** argv){
    string file_name = argv[1];
    size_t lastindex = file_name.find_last_of(".");
    string kf_prefix = file_name.substr(0, lastindex);
    std::string::size_type idx;
    idx = file_name.rfind('.');
    std::string extension = "";

    if(idx != std::string::npos) extension = file_name.substr(idx+1);

    auto * loaded_kf = kDataFrame::load(kf_prefix);

    // Calculating percentile
    auto perc_kf_it = loaded_kf->begin();
    auto* kmerCounts = new(vector<uint64_t>);
    uint64_t real_no_of_kmers = 0;
    uint64_t removed_kmers_from_percentile = 0;

    while (perc_kf_it != loaded_kf->end()) {
        kmerCounts->push_back(perc_kf_it.getCount());
        real_no_of_kmers++;
        perc_kf_it++;
    }
    sort(kmerCounts->begin(), kmerCounts->end());
    uint64_t _idx = (uint64_t)ceil((real_no_of_kmers * LOWEST_PERCENTILE / 100));
    uint64_t count_percentile_cutoff = kmerCounts->at(_idx);
    cout << "   calculated percentile cutoff kmercount=" << count_percentile_cutoff << endl;
    delete kmerCounts;

    kDataFrame * new_kf;
    if(extension == "phmap") new_kf = new kDataFramePHMAP(loaded_kf->getkSize());
    else if(extension == "mqf") new_kf = new kDataFrameMQF(loaded_kf->getkSize());
    else exit (1);

    auto it = loaded_kf->begin();
    while(it != loaded_kf->end()){
        if(it.getCount() <= count_percentile_cutoff){
            removed_kmers_from_percentile++;
            it++;
            continue;
        }
        new_kf->insert(it.getHashedKmer(), it.getCount());
        it++;
    }

    cout << "   Removed kmers from percentile=(" << removed_kmers_from_percentile <<") out of ("<< real_no_of_kmers <<")" << endl;
    new_kf->save(kf_prefix + "_perc" + to_string(LOWEST_PERCENTILE));

}
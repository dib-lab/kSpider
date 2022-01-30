#include <iostream>
#include <cstdint>
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/functional/hash.hpp>

using boost::adaptors::transformed;
using boost::algorithm::join;
using namespace std;
using namespace phmap;

typedef std::chrono::high_resolution_clock Time;

class Combo {

public:
    Combo() = default;

    std::vector<std::pair<uint32_t, uint32_t>> combs;

    void combinations(int n) {
        this->combs.clear();
        this->comb(n, this->r, this->arr);
    }

private:
    int *arr = new int[2];
    int r = 2;

    void comb(int n, int r, int *arr) {
        for (int i = n; i >= r; i--) {
            // choose the first element
            arr[r - 1] = i;
            if (r > 1) { // if still needs to choose
                // recursive into smaller problem
                comb(i - 1, r - 1, arr);

            } else {
                this->combs.emplace_back(std::make_pair(arr[0] - 1, arr[1] - 1));
            }
        }
    }

};

namespace kSpider{
void pairwise(string index_prefix) {

    // Read colors
    string colors_map = index_prefix + "colors.intvectors";
    ifstream input(colors_map.c_str());
    int size;
    input >> size;
    flat_hash_map<uint64_t, vector<uint32_t>> color_to_ids = flat_hash_map<uint64_t, std::vector<uint32_t>>(size);
    for (int i = 0; i < size; i++) {
        uint64_t color, colorSize;
        input >> color >> colorSize;
        uint32_t sampleID;
        color_to_ids[color] = std::vector<uint32_t>(colorSize);
        for (int j = 0; j < colorSize; j++) {
            input >> sampleID;
            color_to_ids[color][j] = sampleID;
        }
    }

    flat_hash_map<uint32_t, uint32_t> colorsCount;

//    auto *ckf = colored_kDataFrame::load(index_prefix);
    auto *kf = kDataFrame::load(index_prefix);
//    auto *kf = ckf->getkDataFrame();

    auto it = kf->begin();
    while (it != kf->end()) {
        colorsCount[it.getCount()]++;
        it++;
    }

    // Free some memory
    delete kf;

    flat_hash_map<uint32_t, uint32_t> groupID_to_kmerCount;
    for (const auto &record : color_to_ids) {
        uint32_t colorCount = colorsCount[record.first];
        for (auto group_id: record.second) {
            groupID_to_kmerCount[group_id] += colorCount;
        }
    }

    std::ofstream fstream_kmerCount;
    fstream_kmerCount.open(index_prefix + "_kSpider_seqToKmersNo.tsv");
    fstream_kmerCount << "ID\tseq\tkmers\n";
    uint64_t counter = 0;
    for(const auto & item : groupID_to_kmerCount){
        fstream_kmerCount << ++counter << '\t' << item.first << '\t' << item.second << '\n';
    }
    fstream_kmerCount.close();

    Combo combo = Combo();
    flat_hash_map<std::pair<uint32_t, uint32_t>, uint32_t, boost::hash<pair<uint32_t, uint32_t>>> edges;
    for (const auto &item : color_to_ids) {
        combo.combinations(item.second.size());
        for (auto const &seq_pair : combo.combs) {
            uint32_t _seq1 = item.second[seq_pair.first];
            uint32_t _seq2 = item.second[seq_pair.second];
            _seq1 > _seq2 ? edges[{_seq1, _seq2}] += colorsCount[item.first] : edges[{_seq2,
                                                                                      _seq1}] += colorsCount[item.first];
        }
    }


    std::ofstream myfile;

    myfile.open(index_prefix + "_kSpider_pairwise.tsv");

    myfile << "ID" << '\t' << "seq1" << '\t' << "seq2" << '\t' << "shared_kmers" << '\n';

    uint64_t line_count = 0;

    for (const auto &edge : edges) {
        myfile << ++line_count << '\t' << edge.first.first << '\t' << edge.first.second << '\t' << edge.second << '\n';
    }

    myfile.close();

}
}
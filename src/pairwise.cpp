#include <iostream>
#include <cstdint>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include "kDataframes/kDataFrameSTL.hpp"
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/functional/hash.hpp>
#include<algorithm>

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
    int* arr = new int[2];
    int r = 2;

    void comb(int n, int r, int* arr) {
        for (int i = n; i >= r; i--) {
            // choose the first element
            arr[r - 1] = i;
            if (r > 1) { // if still needs to choose
                // recursive into smaller problem
                comb(i - 1, r - 1, arr);

            }
            else {
                this->combs.emplace_back(std::make_pair(arr[0] - 1, arr[1] - 1));
            }
        }
    }

};

namespace kSpider {
    void pairwise(string index_prefix) {

        // Read colors
        flat_hash_map<uint32_t, uint32_t> colorsCount;

        auto* kf = kDataFrame::load(index_prefix);
        
        // old_indexing
        // auto* colorColumn = (deduplicatedColumn<StringColorColumn>*) kf->columns["color"];

        // priorityQueue
        auto * colorColumn = (deduplicatedColumn<mixVectors>*) kf->columns["color"];


        flat_hash_map<uint64_t, std::vector<uint32_t>> color_to_ids;

        auto kf_it = kf->begin();
        while (kf_it != kf->end()) {
            uint32_t color_id = colorColumn->index[kf_it.getOrder()];
            colorsCount[color_id]++;
            // old indexing
            // vector<uint32_t> group_ids = colorColumn->values->colors->get(color_id);
            vector<uint32_t> group_ids = kf->getKmerColumnValue<deduplicatedColumn<mixVectors> >("color", kf_it.getKmer());

            color_to_ids[color_id] = std::vector<uint32_t>();
            for (auto& grp_id : group_ids) {
                if (grp_id >= 0) {
                    color_to_ids[color_id].emplace_back(grp_id);
                }
            }
            kf_it++;
        }

        // Free some memory
        cout << "[INFO]" << "deleteing loaded kf." << endl;
        delete kf;


        // Computing unique canonical kmer count from the colors and groups. 
        cout << "[INFO]" << "kmer counting..." << endl;
        flat_hash_map<uint32_t, uint32_t> groupID_to_kmerCount;
        for (auto record : color_to_ids) {
            uint32_t colorCount = colorsCount[record.first];
            for (auto group_id : record.second) {
                groupID_to_kmerCount[group_id] += colorCount;
            }
        }


        std::ofstream fstream_kmerCount;
        fstream_kmerCount.open(index_prefix + "_kSpider_seqToKmersNo.tsv");
        fstream_kmerCount << "ID\tseq\tkmers\n";
        uint64_t counter = 0;
        for (const auto& item : groupID_to_kmerCount) {
            fstream_kmerCount << ++counter << '\t' << item.first << '\t' << item.second << '\n';
        }
        fstream_kmerCount.close();

        Combo combo = Combo();
        flat_hash_map<std::pair<uint32_t, uint32_t>, uint32_t, boost::hash<pair<uint32_t, uint32_t>>> edges;
        for (const auto& item : color_to_ids) {
            combo.combinations(item.second.size());
            for (auto const& seq_pair : combo.combs) {
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

        for (const auto& edge : edges) {
            myfile << ++line_count << '\t' << edge.first.first << '\t' << edge.first.second << '\t' << edge.second << '\n';
        }

        myfile.close();

    }
}
#include <iostream>
#include <cstdint>
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/functional/hash.hpp>
#include <ctime>
#include<omp.h>

using boost::adaptors::transformed;
using boost::algorithm::join;
using namespace std;
using namespace phmap;

using Map = parallel_flat_hash_map<std::pair<uint32_t, uint32_t>, std::uint64_t, boost::hash<pair<uint32_t, uint32_t>>, std::equal_to<std::pair<uint32_t, uint32_t>>, std::allocator<std::pair<const std::pair<uint32_t, uint32_t>, uint32_t>>, 12, std::mutex>;
using int_int_map = parallel_flat_hash_map<uint32_t, uint32_t, std::hash<uint32_t>, std::equal_to<uint32_t>, std::allocator<std::pair<const uint32_t, uint32_t>>, 1>;
using int_vec_map = parallel_flat_hash_map<uint32_t, vector<uint32_t>, std::hash<uint32_t>, std::equal_to<uint32_t>, std::allocator<std::pair<const uint32_t, vector<uint32_t>>>, 1>;

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


template <typename T>
void ascending(T& dFirst, T& dSecond)
{
    if (dFirst > dSecond)
        std::swap(dFirst, dSecond);
}

inline void map_insert(int_int_map& _MAP, uint32_t& key, uint32_t& value) {
    _MAP.insert(make_pair(key, value));
}


namespace kSpider {
    void pairwise(string index_prefix, int user_threads) {

        // Read colors
        auto begin_time = Time::now();
        string colors_map = index_prefix + "colors.intvectors";
        ifstream input(colors_map.c_str());
        int size;
        input >> size;
        int_vec_map color_to_ids;
        for (int i = 0; i < size; i++) {
            uint64_t color, colorSize;
            input >> color >> colorSize;
            uint32_t sampleID;
            color_to_ids.insert(make_pair(color, vector<uint32_t>(colorSize)));
            for (int j = 0; j < colorSize; j++) {
                input >> sampleID;
                color_to_ids[color][j] = sampleID;
            }
        }

        cout << "mapping colors to groups: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;

        begin_time = Time::now();
        int_int_map colorsCount;
        auto* ckf = colored_kDataFrame::load(index_prefix);
        auto namesmap = ckf->names_map();
        auto * kf = ckf->getkDataFrame();
        auto it = kf->begin();
        while (it != kf->end()) {
            colorsCount[it.getCount()]++;
            it++;
        }
        // Free some memory
        delete kf;
        delete ckf;


        // TODO: should be csv, rename later.
        // std::ifstream data(index_prefix + "_kSpider_colorCount.tsv");
        // if (!data.is_open()) std::exit(EXIT_FAILURE);
        // std::string str;
        // std::getline(data, str); // skip the first line
        // while (std::getline(data, str))
        // {
        //     std::istringstream iss(str);
        //     std::string token;
        //     vector<uint32_t> tmp;
        //     while (std::getline(iss, token, ','))
        //         tmp.push_back(stoi(token));
        //     colorsCount.insert(make_pair(tmp[0], tmp[1]));
        // }

        cout << "parsing index colors: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;
        begin_time = Time::now();
        flat_hash_map<uint32_t, uint32_t> groupID_to_kmerCount;
        for (const auto& record : color_to_ids) {
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
        cout << "kmer counting: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;

        begin_time = Time::now();
        clock_t begin_detailed_pairwise_comb, begin_detailed_pairwise_edges, begin_detailed_pairwise_edges_insertion;
        float detailed_pairwise_comb = 0.0;
        float detailed_pairwise_edges = 0.0;
        float detailed_pairwise_edges_insertion = 0.0;

        Map edges;

        // convert map to vec for parallelization purposes.
        auto vec_color_to_ids = std::vector<std::pair<uint32_t, vector<uint32_t>>>(color_to_ids.begin(), color_to_ids.end());

        int thread_num, num_threads, start, end, vec_i;
        int n = vec_color_to_ids.size();

        omp_set_num_threads(user_threads);
        begin_time = Time::now();

#pragma omp parallel private(vec_i,thread_num,num_threads,start,end)
        {
            thread_num = omp_get_thread_num();
            num_threads = omp_get_num_threads();
            start = thread_num * n / num_threads;
            end = (thread_num + 1) * n / num_threads;

            for (vec_i = start; vec_i != end; ++vec_i) {
                auto item = vec_color_to_ids[vec_i];
                Combo combo = Combo();
                combo.combinations(item.second.size());
                for (uint32_t i = 0; i < combo.combs.size(); i++) {
                    // for (auto const& seq_pair : combo.combs) {
                    auto const& seq_pair = combo.combs[i];
                    uint32_t _seq1 = item.second[seq_pair.first];
                    uint32_t _seq2 = item.second[seq_pair.second];
                    ascending(_seq1, _seq2);

                    auto _p = make_pair(_seq1, _seq2);
                    uint32_t ccount = colorsCount[item.first];
                    edges.lazy_emplace_l(_p,
                        [ccount](Map::value_type& v) { v.second++; },           // called only when key was already present
                        [_p, ccount](const Map::constructor& ctor) {
                            ctor(_p, ccount); }
                    ); // construct value_type in place when key not present 
                }
            }
        }

        cout << "pairwise hashmap construction: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;
        cout << "writing pairwise matrix to " << index_prefix << "_kSpider_pairwise.tsv" << endl;

        std::ofstream myfile;
        std::ofstream myfile2;

        myfile.open(index_prefix + "_kSpider_pairwise.tsv");
        myfile2.open(index_prefix + "_kSpider_pairwise_named.tsv");
        myfile << "bin_1" << '\t' << "bin_2" << '\t' << "shared_kmers" << '\t' << "max_containment" << '\n';
        myfile2 << "bin_1" << '\t' << "bin_2" << '\t' << "shared_kmers" << '\t' << "max_containment" << '\n';
        uint64_t line_count = 0;
        for (const auto& edge : edges) {
            float max_containment = (float)edge.second / min(groupID_to_kmerCount[edge.first.first], groupID_to_kmerCount[edge.first.second]);
            myfile << edge.first.first << '\t' << edge.first.second << '\t' << edge.second << '\t' << max_containment << '\n';
            myfile2 << namesmap[edge.first.first] << '\t' << namesmap[edge.first.second] << '\t' << edge.second << '\t' << max_containment << '\n';

        }
        myfile.close();
        myfile2.close();
    }
}
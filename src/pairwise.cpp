#include <iostream>
#include <cstdint>
#include <chrono>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/functional/hash.hpp>
#include <ctime>
#include<omp.h>
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include <cassert>

using boost::adaptors::transformed;
using boost::algorithm::join;
using namespace std;
using namespace phmap;

// using Map = parallel_flat_hash_map<std::pair<uint32_t, uint32_t>, std::uint64_t, boost::hash<pair<uint32_t, uint32_t>>, std::equal_to<std::pair<uint32_t, uint32_t>>, std::allocator<std::pair<const std::pair<uint32_t, uint32_t>, uint32_t>>, 12, std::mutex>;
using int_int_map = parallel_flat_hash_map<uint32_t, uint32_t, std::hash<uint32_t>, std::equal_to<uint32_t>, std::allocator<std::pair<const uint32_t, uint32_t>>, 1>;
using int_vec_map = parallel_flat_hash_map<uint32_t, vector<uint32_t>, std::hash<uint32_t>, std::equal_to<uint32_t>, std::allocator<std::pair<const uint32_t, vector<uint32_t>>>, 1>;

using PAIRS_COUNTER = phmap::parallel_flat_hash_map<
    std::pair<uint32_t, uint32_t>,
    std::uint64_t,
    boost::hash<pair<uint32_t, uint32_t>>,
    std::equal_to<std::pair<uint32_t, uint32_t>>,
    std::allocator<std::pair<const std::pair<uint32_t, uint32_t>, uint64_t>>, 12, std::mutex>;

using BINS_KMER_COUNT = phmap::parallel_flat_hash_map<
    std::string, uint32_t,
    phmap::priv::hash_default_hash<std::string>,
    phmap::priv::hash_default_eq<std::string>,
    std::allocator<std::pair<const std::string, uint32_t>>,
    1,
    std::mutex>;

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

    void set_to_vector(const phmap::flat_hash_set<uint32_t>& set, vector<uint32_t>& vec) {
        vec.clear();
        vec.reserve(set.size());
        for (auto& i : set) {
            vec.push_back(i);
        }
    }

    void load_colors_to_sources(const std::string& filename, int_vec_map * map)
    {
        phmap::BinaryInputArchive ar_in(filename.c_str());
        size_t size;
        ar_in.loadBinary(&size);
        map->reserve(size);
        while (size--)
        {
            uint64_t k;
            phmap::flat_hash_set<uint32_t> v;
            vector<uint32_t> vVec;
            set_to_vector(v, vVec);
            ar_in.loadBinary(&k);
            ar_in.loadBinary(&v);
            map->insert_or_assign(std::move(k), std::move(vVec));
        }
    }

    void load_colors_count(const std::string& filename, int_int_map& map) {
        flat_hash_map<uint64_t, uint64_t> tmpMap;
        phmap::BinaryInputArchive ar_in_colorsCount(filename.c_str());
        tmpMap.phmap_load(ar_in_colorsCount);
        assert(tmpMap.size());
        for (auto& i : tmpMap) {
            map.insert_or_assign(i.first, i.second);
        }
    }

    void pairwise(string index_prefix, int user_threads) {

        // Read colors

        int_vec_map color_to_ids; // = new phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint32_t>>;
        string colors_map_file = index_prefix + "_color_to_sources.bin";
        load_colors_to_sources(colors_map_file, &color_to_ids);

        auto begin_time = Time::now();

        cout << "mapping colors to groups: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;

        begin_time = Time::now();
        int_int_map colorsCount;

        load_colors_count(index_prefix + "_color_count.bin", colorsCount);

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

        PAIRS_COUNTER edges;

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
                        [ccount](PAIRS_COUNTER::value_type& v) { v.second++; },           // called only when key was already present
                        [_p, ccount](const PAIRS_COUNTER::constructor& ctor) {
                            ctor(_p, ccount); }
                    ); // construct value_type in place when key not present 
                }
            }
        }

        cout << "pairwise hashmap construction: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;
        cout << "writing pairwise matrix to " << index_prefix << "_kSpider_pairwise.tsv" << endl;

        std::ofstream myfile;
        myfile.open(index_prefix + "_kSpider_pairwise.tsv");
        myfile
            << "source_1"
            << "\tsource_2"
            << "\tshared_kmers"
            << "\tmin_containment"
            << "\tavg_containment"
            << "\tmax_containment"
            << '\n';
        uint64_t line_count = 0;
        for (const auto& edge : edges) {
            uint64_t shared_kmers = edge.second;
            uint32_t source_1 = edge.first.first;
            uint32_t source_2 = edge.first.second;
            uint32_t source_1_kmers = groupID_to_kmerCount[source_1];
            uint32_t source_2_kmers = groupID_to_kmerCount[source_2];

            float cont_1_in_2 = (float)shared_kmers / source_2_kmers;
            float cont_2_in_1 = (float)shared_kmers / source_1_kmers;
            float min_containment = min(cont_1_in_2, cont_2_in_1);
            float avg_containment = (cont_1_in_2 + cont_2_in_1) / 2.0;
            float max_containment = max(cont_1_in_2, cont_2_in_1);

            myfile
                << source_1
                << '\t' << source_2
                << '\t' << shared_kmers
                << '\t' << min_containment
                << '\t' << avg_containment
                << '\t' << max_containment
                << '\n';
        }
        myfile.close();
    }
}
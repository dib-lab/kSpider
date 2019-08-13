#include "kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "stdint.h"
#include <vector>
#include <string>
#include <set>
#include "pairwise_matrix.hpp"
#include <boost/functional/hash.hpp>

using phmap::flat_hash_map;
using std::vector;
using std::set;
using std::string;


// A hash function used to hash a pair of any kind
struct hash_pair {
    template<class T1, class T2>
    size_t operator()(const pair<T1, T2> &p) const {
        return std::hash<T1>{}(p.first) ^ std::hash<T2>{}(p.second);
    }
};


class virtualQs {
public:
    flat_hash_map<uint8_t, set<uint64_t> > temp_superColors;
    flat_hash_map<uint8_t, flat_hash_map<uint64_t, set<uint64_t> >> superColors;
    flat_hash_map<uint8_t, flat_hash_map<uint64_t, uint64_t> > superColorsCount;
    flat_hash_map<int, uint64_t> masks;
//    vector<flat_hash_map<pair<uint32_t ,uint32_t>, uint32_t, boost::hash<pair<uint32_t, uint32_t>> >> edges;
//    flat_hash_map<std::pair<uint32_t, uint32_t>, flat_hash_map<uint8_t , uint32_t>, hash_pair> edges;
    flat_hash_map<std::pair<pair<uint32_t, uint32_t>, uint8_t>, uint32_t, boost::hash<pair<pair<uint32_t, uint32_t>, uint8_t>>> edges;
//    flat_hash_map<uint32_t, flat_hash_map<uint32_t, flat_hash_map<uint8_t, uint16_t>>> edges;
    flat_hash_map<uint32_t, uint32_t> seq_to_kmers_no;
    flat_hash_map<uint64_t, vector<uint32_t>> color_to_ids;
//    flat_hash_map<uint32_t, uint32_t> namesMap;

    vector<int> mainQs;

    kDataFrame *KF;
    string index_prefix;
    int kSize;

    virtualQs(string index_path, uint8_t minQ, uint8_t maxQ, uint8_t stepQ);

    uint64_t create_super_color(set<uint64_t> &colors);

    void calculate_kmers_number();

    void pairwise();


};
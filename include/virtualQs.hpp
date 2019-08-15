#include "kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "stdint.h"
#include <vector>
#include <string>
#include <set>
#include <sqlite3.h>
#include <boost/functional/hash.hpp>

using phmap::flat_hash_map;
using phmap::flat_hash_set;
using std::vector;
using std::set;
using std::string;


class virtualQs {
public:
    flat_hash_map<uint8_t, flat_hash_set<uint64_t> > temp_superColors;
    flat_hash_map<uint8_t, flat_hash_map<uint64_t, flat_hash_set<uint64_t> >> superColors;
    flat_hash_map<uint8_t, flat_hash_map<uint64_t, uint64_t> > superColorsCount;
    flat_hash_map<int, uint64_t> masks;
//    vector<flat_hash_map<pair<uint32_t ,uint32_t>, uint32_t, boost::hash<pair<uint32_t, uint32_t>> >> edges;


    //1
//    flat_hash_map<std::pair<uint32_t, uint32_t>, flat_hash_map<uint8_t, uint32_t>, boost::hash<pair<uint32_t, uint32_t>>> edges;

    //2
        flat_hash_map<std::pair<pair<uint32_t, uint32_t>, uint8_t>, uint32_t, boost::hash<pair<pair<uint32_t, uint32_t>, uint8_t>>> edges;




//    flat_hash_map<uint32_t, flat_hash_map<uint32_t, flat_hash_map<uint8_t, uint16_t>>> edges;
    flat_hash_map<uint32_t, uint32_t> seq_to_kmers_no;
    flat_hash_map<uint64_t, vector<uint32_t>> color_to_ids;
//    flat_hash_map<uint32_t, uint32_t> namesMap;

    sqlite3 *DB;
    char *DB_ErrMsg = 0;

    set<int > mainQs;

    kDataFrame *KF;
    string index_prefix;
    int kSize;

    virtualQs(string index_path, set<int> mainQs);

    uint64_t create_super_color(flat_hash_set<uint64_t> &colors);

    void calculate_kmers_number();

    void pairwise();

    void export_to_sqlite();
    void export_to_tsv();

};
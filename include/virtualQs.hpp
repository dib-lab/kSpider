#include "kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "stdint.h"
#include <vector>
#include <string>
#include <set>

using phmap::flat_hash_map;
using std::vector;
using std::set;
using std::string;

class virtualQs {
public:
    flat_hash_map<uint8_t, set<uint64_t> > temp_superColors;
    flat_hash_map<uint8_t, flat_hash_map<uint64_t, set<uint64_t> >> superColors;
    flat_hash_map<uint8_t, flat_hash_map<uint64_t, uint64_t> > superColorsCount;
    flat_hash_map<int, uint64_t> masks;
    vector<int> mainQs;

    kDataFrame *KF;
    string index_prefix;
    int kSize;

    virtualQs(string index_path, uint8_t minQ, uint8_t maxQ, uint8_t stepQ);

    uint64_t create_super_color(set<uint64_t > &colors);


private:

    void prepare();


};
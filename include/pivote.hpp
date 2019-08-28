#include "parallel_hashmap/phmap.h"
#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include <stdint.h>
#include <sstream>
#include <boost/functional/hash.hpp>
#include <set>

using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::string;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::set;
using phmap::flat_hash_map;

typedef flat_hash_map<pair<pair<uint32_t, uint32_t>, uint16_t>, vector<pair<uint8_t, uint16_t>>, boost::hash<pair<pair<uint32_t, uint32_t>, uint16_t>>> MAP;
typedef flat_hash_map<pair<pair<uint32_t, uint32_t>, uint16_t>, flat_hash_map<uint8_t, uint16_t>, boost::hash<pair<pair<uint32_t, uint32_t>, uint16_t>>> MAP2;

void pivote(string index_prefex, set<int> allQs);
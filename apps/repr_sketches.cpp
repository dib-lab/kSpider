#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <cstdint>
#include <unordered_map>
#include <parallel_hashmap/phmap.h>

using namespace boost::algorithm;
using namespace std;

bool comp(pair<uint64_t,uint64_t> a, pair<uint64_t,uint64_t> b) {
    return a.second > b.second;
}


int main(int argc, char** argv) {
    ifstream fin(argv[1]);
    phmap::flat_hash_map<uint64_t, uint64_t> count;
    string line;
    getline(fin, line); // skip header.
    while (getline(fin, line)) {
        // Split line into tab-separated parts
        vector<string> parts;
        split(parts, line, boost::is_any_of("\t"));
        float containment = stof(parts[4]);
        if (containment > 0.20) {
            uint64_t from_node = stoi(parts[0]);
            uint64_t to_node = stoi(parts[1]);
            count[from_node]++;
            count[to_node]++;
        }

    }
    fin.close();

    std::vector<std::pair<uint64_t, uint64_t>> elems(count.begin(), count.end());
    std::sort(elems.begin(), elems.end(), comp);

    for (auto& [k, v] : elems) {
        cout << k << ": " << v << endl;
    }
}
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

int main(int argc, char** argv) {
    string index_prefix = argv[1];

    flat_hash_map<uint32_t, uint32_t> colorsCount;
    auto* kf = kDataFrame::load(index_prefix);
    auto it = kf->begin();
    while (it != kf->end()) {
        uint32_t ccount = it.getCount();
        colorsCount[ccount]++;
        it++;
    }

    std::ofstream fstream_ccount;
    fstream_ccount.open(index_prefix + "_kSpider_colorCount.tsv");
    fstream_ccount << "color,count\n";
    for (auto [color, count] : colorsCount)
        fstream_ccount << color << "," << count << "\n";
    fstream_ccount.close();

}
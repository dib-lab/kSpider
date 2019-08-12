#include "virtualQs.hpp"
#include <iostream>
#include <bitset>

using std::cerr;
using std::endl;

inline uint64_t create_mask(unsigned kSize, unsigned Q) {
    return ((1ULL << Q * 2ULL) - 1ULL) << (kSize * 2ULL - Q * 2ULL);
}

virtualQs::virtualQs(string index_prefix, uint8_t minQ, uint8_t maxQ, uint8_t stepQ) {

    // Loading the index
    cerr << "[INFO] Loading the index" << endl;
    this->KF = kDataFrame::load(index_prefix);
    this->kSize = KF->ksize();
    this->index_prefix = index_prefix;

    // Constructing masks
    for (int Q = minQ; Q <= maxQ; Q += stepQ)
        this->mainQs.push_back(Q);

    bool ksize_in_Qs = (this->masks.find(this->kSize) != this->masks.end());

    if (!ksize_in_Qs)
        this->mainQs.push_back(this->kSize);

    for (auto const &Q : this->mainQs) {
        this->masks[Q] = create_mask(this->kSize, Q);
        this->superColors[Q] = flat_hash_map<uint64_t, set<uint64_t>>();
        this->superColorsCount[Q] = flat_hash_map<uint64_t, uint64_t>();
        this->temp_superColors[Q] = set<uint64_t>();
    }
}

uint64_t virtualQs::create_super_color(set<uint64_t> &colors) {
    uint64_t seed = colors.size();
    for (auto &color : colors)
        seed ^= color + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}
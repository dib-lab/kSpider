#include "virtualQs.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include "combinations.hpp"

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

    string colors_map = this->index_prefix + "colors.intvectors";
    ifstream input(colors_map.c_str());
    int size;
    input >> size;
    color_to_ids = flat_hash_map<uint64_t, std::vector<uint32_t>>(size);
    for (int i = 0; i < size; i++) {
        uint64_t color, colorSize;
        input >> color >> colorSize;
        uint32_t sampleID;
        color_to_ids[color] = std::vector<uint32_t>(colorSize);
        for (int j = 0; j < colorSize; j++) {
            input >> sampleID;
            color_to_ids[color][j] = sampleID;
        }
    }


    // Read NamesMap

//    ifstream namesMapIn(index_prefix + ".namesMap");
//    namesMapIn >> size;
//    for (int i = 0; i < size; i++) {
//        uint32_t sample_id;
//        string sample_name;
//        namesMapIn >> sample_id >> sample_name;
//        this->namesMap[sample_id] = sample_id;
//    }

}

uint64_t virtualQs::create_super_color(set<uint64_t> &colors) {
    uint64_t seed = colors.size();
    for (auto &color : colors)
        seed ^= color + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

void virtualQs::calculate_kmers_number() {

    for (auto const &superColor : this->superColors[this->kSize]) {
        vector<uint32_t> tr_ids;
        for (auto const &color : superColor.second)
            for (auto const &id : this->color_to_ids[color])
                tr_ids.push_back(id);

        uint32_t color_count = this->superColorsCount[this->kSize][superColor.first];

        for (auto const &tr_id : tr_ids) {
            bool tr_exist = (this->seq_to_kmers_no.find(tr_id) != this->seq_to_kmers_no.end());
            if (tr_exist) {
                this->seq_to_kmers_no[tr_id] += color_count;
            } else {
                this->seq_to_kmers_no[tr_id] = color_count;
            }
        }

    }

}

void virtualQs::pairwise() {
    this->calculate_kmers_number();
    Combo combo = Combo();
    for (auto const &Q : this->mainQs) {
        cerr << "Processing Q: " << Q << endl;
        for (auto const &superColor : this->superColors[Q]) {
            vector<uint32_t> tr_ids;
            for (auto const &color : superColor.second)
                for (auto const &id : this->color_to_ids[color])
                    tr_ids.push_back(id);

            uint32_t color_count = this->superColorsCount[Q][superColor.first];
            combo.combinations(tr_ids.size());
            for (auto const &seq_pair : combo.combs) {
                uint32_t _seq1 = tr_ids[seq_pair.first];
                uint32_t _seq2 = tr_ids[seq_pair.second];
                _seq1 > _seq2 ? this->edges[{{_seq1, _seq2}, Q}] += color_count : this->edges[{{_seq2, _seq1},
                                                                                               Q}] += color_count;
            }
        }

    }
}


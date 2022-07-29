#include <iostream>
#include <cstdint>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include <ctime>
#include<omp.h>
#include "RSJparser.tcc"
#include <glob.h>
#include <string>
#include <stdexcept>
#include "parallel_hashmap/phmap_dump.h"
#include <cstdlib>

using namespace std;
// using namespace phmap;
using JSON = RSJresource;

typedef std::chrono::high_resolution_clock Time;

inline uint64_t unrolled(std::string const& value) {
    uint64_t result = 0;

    size_t const length = value.size();
    switch (length) {
    case 20:    result += (value[length - 20] - '0') * 10000000000000000000ULL;
    case 19:    result += (value[length - 19] - '0') * 1000000000000000000ULL;
    case 18:    result += (value[length - 18] - '0') * 100000000000000000ULL;
    case 17:    result += (value[length - 17] - '0') * 10000000000000000ULL;
    case 16:    result += (value[length - 16] - '0') * 1000000000000000ULL;
    case 15:    result += (value[length - 15] - '0') * 100000000000000ULL;
    case 14:    result += (value[length - 14] - '0') * 10000000000000ULL;
    case 13:    result += (value[length - 13] - '0') * 1000000000000ULL;
    case 12:    result += (value[length - 12] - '0') * 100000000000ULL;
    case 11:    result += (value[length - 11] - '0') * 10000000000ULL;
    case 10:    result += (value[length - 10] - '0') * 1000000000ULL;
    case  9:    result += (value[length - 9] - '0') * 100000000ULL;
    case  8:    result += (value[length - 8] - '0') * 10000000ULL;
    case  7:    result += (value[length - 7] - '0') * 1000000ULL;
    case  6:    result += (value[length - 6] - '0') * 100000ULL;
    case  5:    result += (value[length - 5] - '0') * 10000ULL;
    case  4:    result += (value[length - 4] - '0') * 1000ULL;
    case  3:    result += (value[length - 3] - '0') * 100ULL;
    case  2:    result += (value[length - 2] - '0') * 10ULL;
    case  1:    result += (value[length - 1] - '0');
    }
    return result;
}

template<>
uint64_t RSJresource::as<uint64_t>(const uint64_t& def) {
    if (!exists()) return (0); return (unrolled(data));
}


int main(int argc, char** argv) {

    if (argc != 5) {
        cout << "run: ./sig_to_bin <sig> <kSize> <min_abundance> <output_file>" << endl;
        exit(1);
    }

    string sig_path = argv[1];
    int kSize = stoi(argv[2]);
    int min_abundance = stoi(argv[3]);
    string output_path = argv[4];

    auto begin_time = Time::now();
    zstr::ifstream sig_stream(sig_path);
    JSON sig(sig_stream);
    phmap::flat_hash_set<uint64_t> tmp_hashes;
    int number_of_sub_sigs = sig[0]["signatures"].size();
    for (int i = 0; i < number_of_sub_sigs; i++) {
        int current_kSize = sig[0]["signatures"][i]["ksize"].as<int>();
        auto loaded_sig_it = sig[0]["signatures"][i]["mins"].as_array().begin();
        auto abundance_it = sig[0]["signatures"][i]["abundances"].as_array().begin();
        if (current_kSize == kSize) {
            while (loaded_sig_it != sig[0]["signatures"][i]["mins"].as_array().end()) {
                int abund = abundance_it->as<int>();
                if (abund >= min_abundance)
                    tmp_hashes.insert(loaded_sig_it->as<uint64_t>());
                loaded_sig_it++;
                abundance_it++;
            }
            break;
        }
    }


    cout << "inserted " << tmp_hashes.size() << " hashes." << endl;
    string out_path = output_path;
    phmap::BinaryOutputArchive ar_out(out_path.c_str());
    tmp_hashes.phmap_dump(ar_out);
    cout << "Conversion done in " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;


}
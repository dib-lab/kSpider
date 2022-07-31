#include <iostream>
#include <cstdint>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include <ctime>
#include <omp.h>
#include "cpp-json/json.h"
#include "zstr.hpp"
#include <glob.h>
#include <string>
#include <stdexcept>
#include "parallel_hashmap/phmap_dump.h"
#include <cstdlib>

using namespace std;
// using namespace phmap;

typedef std::chrono::high_resolution_clock Time;


int main(int argc, char** argv) {

    if (argc != 3) {
        cout << "run: ./dump_sig <sig> <kSize>" << endl;
        exit(1);
    }

    string sig_path = argv[1];
    int kSize = stoi(argv[2]);

    phmap::flat_hash_set<uint64_t> tmp_hashes;

    auto begin_time = Time::now();
    zstr::ifstream sig_stream(sig_path);
    json::value json = json::parse(sig_stream);
    auto sourmash_sig = json[0]["signatures"];
    const json::array& sig_array = as_array(sourmash_sig);
    for (auto it = sig_array.begin(); it != sig_array.end(); ++it) {
        const json::value& v = *it;
        if (v["ksize"] == kSize) {
            const json::array& mins = as_array(v["mins"]);
            auto mins_it = mins.begin();
            while (mins_it != mins.end()) {
                tmp_hashes.insert(json::to_number<uint64_t>(*mins_it));
                mins_it++;
            }
        }
        break;
    }


    for (const uint64_t& hash : tmp_hashes) cout << hash << endl;


}
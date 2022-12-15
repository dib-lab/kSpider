#include <iostream>
#include <cstdint>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include <ctime>
#include <omp.h>
#include <glob.h>
#include <string>
#include <stdexcept>
#include "parallel_hashmap/phmap_dump.h"
#include <cstdlib>
#include "cpp-json/json.h"
#include "zstr.hpp"

using namespace std;
// using namespace phmap;

typedef std::chrono::high_resolution_clock Time;


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

    phmap::flat_hash_set<uint64_t> tmp_hashes;

    zstr::ifstream sig_stream(sig_path);
    json::value json = json::parse(sig_stream);
    auto sourmash_sig = json[0]["signatures"];
    const json::array& sig_array = as_array(sourmash_sig);
    for (auto it = sig_array.begin(); it != sig_array.end(); ++it) {
        const json::value& v = *it;
        if (v["ksize"] == kSize) {
            const json::array& mins = as_array(v["mins"]);
            const json::array& abundances = as_array(v["abundances"]);
            auto mins_it = mins.begin();
            auto abund_it = abundances.begin();
            while (mins_it != mins.end()) {
                if (json::to_number<int>(*abund_it) >= min_abundance)
                    tmp_hashes.insert(json::to_number<uint64_t>(*mins_it));

                mins_it++;
                abund_it++;
            }
        }
        break;
    }


    cout << "inserted " << tmp_hashes.size() << " hashes." << endl;
    string out_path = output_path;
    phmap::BinaryOutputArchive ar_out(out_path.c_str());
    tmp_hashes.phmap_dump(ar_out);
    cout << "Conversion done in " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;


}
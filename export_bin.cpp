#include <iostream>
#include <cstdint>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include <ctime>
#include<omp.h>
#include <glob.h>
#include <string>
#include <stdexcept>
#include "parallel_hashmap/phmap_dump.h"
#include <cstdlib>

using namespace std;
// using namespace phmap;


int main(int argc, char** argv) {

    if (argc != 2) {
        cout << "run: ./dump_bin <bin>" << endl;
        exit(1);
    }
    
    string bin_path = argv[1];

    phmap::flat_hash_set<uint64_t> table_in;
    phmap::BinaryInputArchive ar_in(bin_path.c_str());
    table_in.phmap_load(ar_in);
    cerr << "loaded bin size: " << table_in.size() << endl;

    for(const uint64_t & hash : table_in) cout << hash << endl;

}
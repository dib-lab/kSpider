#include <iostream>
#include <kDataFrame.hpp>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithms.hpp>
#include "assert.h"
#include <cstring>
#include <cstdint>
#include <math.h>
#include <fstream>

using namespace std;

int main(int argc, char** argv) {

    int kSize = 25;
    auto* kf = new kDataFramePHMAP(KMERS, mumur_hasher, { {"kSize", kSize} });

    string line;
    string kmc_file = argv[1];

    std::string base_filename = kmc_file.substr(kmc_file.find_last_of("/\\") + 1);
    base_filename = base_filename.substr(0, base_filename.find('_'));

    ifstream FileHandler(kmc_file.c_str());
    while (std::getline(FileHandler, line)) {
        std::vector<string> tokens;
        std::istringstream iss(line);
        std::string token;


        while (std::getline(iss, token, '\t'))
            tokens.push_back(token);


        string kmer = tokens[0];
        uint64_t count;
        std::istringstream _iss(tokens[1]);
        _iss >> count;
        kf->insert(kmer, count);
    }
    kf->save(base_filename);
}
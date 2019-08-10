#include "kCluster.hpp"
#include <iostream>
#include "stdint.h"

int main(int argc, char ** argv){

    std::cout << "Hello World" << std::endl;

    kDataFrame *KF = kDataFrame::load(argv[1]);

    auto it = KF->begin();
    uint64_t kmer;
    uint64_t color;

    while (it != KF->end()){
        kmer = it.getHashedKmer();
        color = it.getKmerCount();
        it++;
    }

    std::cout << KF->ksize() << std::endl;

    return 0;
}
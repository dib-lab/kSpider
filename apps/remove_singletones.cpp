#include "kDataFrame.hpp"
#include "algorithms.hpp"
#include <iostream>
#include <fstream>
using namespace std;

int main(int args, char ** argv){
    string file_name = argv[1];
    size_t lastindex = file_name.find_last_of(".");
    string kf_prefix = file_name.substr(0, lastindex);
    std::string::size_type idx;
    idx = file_name.rfind('.');
    std::string extension = "";
    int removed_singletones = 0;

    if(idx != std::string::npos) extension = file_name.substr(idx+1);

    auto * loaded_kf = kDataFrame::load(kf_prefix);

    kDataFrame * new_kf;
    if(extension == "phmap") new_kf = new kDataFramePHMAP(loaded_kf->getkSize());
    else if(extension == "mqf") new_kf = new kDataFrameMQF(loaded_kf->getkSize());

    uint64_t real_no_of_kmers = 0;

    auto it = loaded_kf->begin();
    while(it != loaded_kf->end()){
        real_no_of_kmers++;
        if(it.getCount() == 1){
            removed_singletones++;
            it++;
            continue;
        }
        new_kf->insert(it.getHashedKmer(), it.getCount());
        it++;
    }

    cout << "   Removed singletones=(" << removed_singletones <<") out of ("<< real_no_of_kmers <<")" << endl;
    new_kf->save(kf_prefix + ".singlefree");

}
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
    if(idx != std::string::npos) extension = file_name.substr(idx+1);

    auto * kf = kDataFrame::load(kf_prefix);
    auto it = kf->begin();

    while(it != kf->end()){
        cout << it.getHashedKmer() << endl;
        it++;
    }

}
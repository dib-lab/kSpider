#include "kDataFrame.hpp"
#include "algorithms.hpp"
#include <iostream>
#include <fstream>
using namespace std;


int main(int args, char** argv) {
    string file_name = argv[1];
    size_t lastindex = file_name.find_last_of(".");
    string kf_prefix = file_name.substr(0, lastindex);
    std::string::size_type idx;
    idx = file_name.rfind('.');
    std::string extension = "";


    auto* loaded_kf = kDataFrame::load(kf_prefix);

    ofstream COUNT(kf_prefix + ".count");

    auto it = loaded_kf->begin();
    while (it != loaded_kf->end()) {
        COUNT << it.getCount() << endl;
        it++;
    }

  // Close the file
  COUNT.close();

}
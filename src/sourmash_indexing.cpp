#include "kSpider.hpp"
#include <iostream>
#include <cstdint>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include "algorithms.hpp"
#include <glob.h>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <string.h>
#include "defaultColumn.hpp"
#include "kDataframes/kDataFrameSTL.hpp"
#include <fstream>

#include "RSJparser.tcc.hpp"
#include <fstream>
#include "zstr.hpp"

using JSON = RSJresource;

// thanks to http://jsteemann.github.io/blog/2016/06/02/fastest-string-to-uint64-conversion-method/
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
    if (!exists()) return (0); // required
    return (unrolled(data)); // example
}

// thanks to https://stackoverflow.com/a/8615450/3371177
inline std::vector<std::string> glob2(const std::string& pattern) {
    using namespace std;

    // glob struct resides on the stack
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // do the glob operation
    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (return_value != 0) {
        globfree(&glob_result);
        stringstream ss;
        ss << "glob() failed with return_value " << return_value << endl;
        throw std::runtime_error(ss.str());
    }

    // collect all the filenames into a std::list<std::string>
    vector<string> filenames;
    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
        filenames.push_back(string(glob_result.gl_pathv[i]));
    }

    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}


namespace kSpider {

    void sourmash_index_kp1(string sigs_dir, int selective_kSize) {

        kDataFrame* frame;
        std::string dir_prefix = sigs_dir.substr(sigs_dir.find_last_of("/\\") + 1);

        flat_hash_map<string, string> namesMap;
        string names_fileName = sigs_dir;

        flat_hash_map<string, uint64_t> tagsMap;
        flat_hash_map<string, uint64_t> groupNameMap;
        auto* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
        flat_hash_map<uint64_t, uint64_t> colorsCount;
        uint64_t readID = 0, groupID = 1;
        string seqName, groupName;
        string line;
        priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
        flat_hash_map<string, uint64_t> groupCounter;

        int total_sigs_number = 0;
        int kframe_kSize = selective_kSize;
        if (selective_kSize > 31) kframe_kSize = 31;
        frame = new kDataFramePHMAP(kframe_kSize, mumur_hasher);



        auto* colors = new deduplicatedColumn< StringColorColumn>();
        frame->addColumn("color", (Column*)colors);

        for (const auto& dirEntry : glob2(sigs_dir + "/*")) {
            string file_name = (string)dirEntry;
            size_t lastindex = file_name.find_last_of(".");
            string sig_prefix = file_name.substr(0, lastindex);


            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if (idx != std::string::npos) extension = file_name.substr(idx + 1);
            if (extension != "sig" && extension != "gz") continue;

            zstr::ifstream tmp_stream(file_name);
            JSON sig(tmp_stream);
            int number_of_sub_sigs = sig[0]["signatures"].size();
            string general_name = sig[0]["name"].as<std::string>();
            if (general_name == "") {
                std::string sig_basename = sig_prefix.substr(sig_prefix.find_last_of("/\\") + 1);
                general_name = sig_basename;
            }

            // uncomment if if groupname = sig_name
            // total_sigs_number++;

            for (int i = 0; i < number_of_sub_sigs; i++) {
                int current_kSize = sig[0]["signatures"][i]["ksize"].as<int>();
                if (current_kSize != selective_kSize) continue;

                total_sigs_number++;


                // std::string sig_basename = sig_prefix.substr(sig_prefix.find_last_of("/\\") + 1);
                string md5sum = sig[0]["signatures"][i]["md5sum"].as<std::string>();
                string sig_name = md5sum + ":" + general_name;

                // Here we can decide
                seqName = sig_name;
                groupName = general_name;

                namesMap.insert(make_pair(seqName, groupName));
                auto it = groupNameMap.find(groupName);
                groupCounter[groupName]++;
                if (it == groupNameMap.end()) {
                    groupNameMap.insert(make_pair(groupName, groupID));
                    tagsMap.insert(make_pair(to_string(groupID), groupID));
                    vector<uint32_t> tmp;
                    tmp.clear();
                    tmp.push_back(groupID);
                    legend->insert(make_pair(groupID, tmp));
                    colorsCount.insert(make_pair(groupID, 0));
                    groupID++;
                }
            }
        }

        cout << "namesmap construction done..." << endl;


        flat_hash_map<uint64_t, string> inv_groupNameMap;
        for (auto& _ : groupNameMap)
            inv_groupNameMap[_.second] = _.first;

        string kmer;
        int __batch_count = 0;
        readID = 0;
        int processed_sigs_count = 0;

        for (const auto& dirEntry : glob2(sigs_dir + "/*")) {
            string file_name = (string)dirEntry;
            size_t lastindex = file_name.find_last_of(".");
            string sig_prefix = file_name.substr(0, lastindex);

            std::string sig_basename = sig_prefix.substr(sig_prefix.find_last_of("/\\") + 1);


            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if (idx != std::string::npos) extension = file_name.substr(idx + 1);
            if (extension != "sig" && extension != "gz") continue;


            zstr::ifstream sig_stream(file_name);
            JSON sig(sig_stream);
            int number_of_sub_sigs = sig[0]["signatures"].size();
            string general_name = sig[0]["name"].as<std::string>();
            if (general_name == "") {
                std::string sig_basename = sig_prefix.substr(sig_prefix.find_last_of("/\\") + 1);
                general_name = sig_basename;
            }

            // start

            for (int i = 0; i < number_of_sub_sigs; i++) {

                int current_kSize = sig[0]["signatures"][i]["ksize"].as<int>();
                if (current_kSize != selective_kSize) continue;

                cout << "Processing " << ++processed_sigs_count << "/" << total_sigs_number << " | " << general_name << " k:" << selective_kSize << " ... " << endl;

                string md5sum = sig[0]["signatures"][i]["md5sum"].as<std::string>();
                string sig_name = md5sum + ":" + general_name;
                string readName = sig_name;
                string groupName = general_name;

                flat_hash_map<uint64_t, uint64_t> convertMap;
                auto loaded_sig_it = sig[0]["signatures"][i]["mins"].as_array().begin();
                uint64_t readTag = groupNameMap.find(groupName)->second;
                convertMap.clear();
                convertMap.insert(make_pair(0, readTag));
                convertMap.insert(make_pair(readTag, readTag));

                while (loaded_sig_it != sig[0]["signatures"][i]["mins"].as_array().end()) {
                    string groupName = sig_basename;
                    uint64_t hashed_kmer = loaded_sig_it->as<uint64_t>();

                    frame->insert(hashed_kmer);

                    uint64_t kmerOrder = frame->getkmerOrder(hashed_kmer);
                    if (colors->size() < frame->size() + 1)
                        colors->resize(frame->size() + 1);
                    uint64_t currentTag = colors->index[kmerOrder];
                    auto itc = convertMap.find(currentTag);
                    if (itc == convertMap.end()) {
                        vector<uint32_t> colors = legend->find(currentTag)->second;
                        auto tmpiT = find(colors.begin(), colors.end(), readTag);
                        if (tmpiT == colors.end()) {
                            colors.push_back(readTag);
                            sort(colors.begin(), colors.end());
                        }

                        string colorsString = to_string(colors[0]);
                        for (unsigned int k = 1; k < colors.size(); k++) {
                            colorsString += ";" + to_string(colors[k]);
                        }

                        auto itTag = tagsMap.find(colorsString);
                        if (itTag == tagsMap.end()) {
                            uint64_t newColor;
                            if (freeColors.empty()) {
                                newColor = groupID++;
                            }
                            else {
                                newColor = freeColors.top();
                                freeColors.pop();
                            }

                            tagsMap.insert(make_pair(colorsString, newColor));
                            legend->insert(make_pair(newColor, colors));
                            itTag = tagsMap.find(colorsString);
                            colorsCount[newColor] = 0;
                        }
                        uint64_t newColor = itTag->second;

                        convertMap.insert(make_pair(currentTag, newColor));
                        itc = convertMap.find(currentTag);
                    }

                    if (itc->second != currentTag) {

                        colorsCount[currentTag]--;
                        if (colorsCount[currentTag] == 0 && currentTag != 0) {
                            auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                            if (_invGroupNameIT == inv_groupNameMap.end()) {
                                freeColors.push(currentTag);

                                vector<uint32_t> colors = legend->find(currentTag)->second;
                                string colorsString = to_string(colors[0]);
                                for (unsigned int k = 1; k < colors.size(); k++) {
                                    colorsString += ";" + to_string(colors[k]);
                                }
                                tagsMap.erase(colorsString);

                                legend->erase(currentTag);
                                if (convertMap.find(currentTag) != convertMap.end())
                                    convertMap.erase(currentTag);
                            }

                        }
                        colorsCount[itc->second]++;
                    }

                    colors->index[kmerOrder] = itc->second;


                    readID += 1;
                    groupCounter[groupName]--;
                    if (colorsCount[readTag] == 0) {
                        if (groupCounter[groupName] == 0) {
                            freeColors.push(readTag);
                            legend->erase(readTag);
                        }
                    }
                    loaded_sig_it++;
                }

            }
        }

        colors->values = new StringColorColumn(legend, groupCounter.size());
        delete legend;
        for (auto& iit : namesMap) {
            uint32_t sampleID = groupNameMap[iit.second];
            colors->values->namesMap[sampleID] = iit.second;
        }

        cout << "saving to " << dir_prefix << " ..." << endl;
        frame->save(dir_prefix);

        ofstream f_namesmap;
        f_namesmap.open(dir_prefix + ".kSpider_namesMap");
        for (std::pair<const unsigned int, class std::basic_string<char> >& name : colors->values->namesMap)
            f_namesmap << name.first << "," << name.second << endl;
        f_namesmap.close();

    }

    void sourmash_index_kp1_fast(string sigs_dir, int selective_kSize) {

        kDataFrame* frame;
        std::string dir_prefix = sigs_dir.substr(sigs_dir.find_last_of("/\\") + 1);

        flat_hash_map<string, string> namesMap;
        string names_fileName = sigs_dir;

        flat_hash_map<string, uint64_t> tagsMap;
        flat_hash_map<string, uint64_t> groupNameMap;
        auto* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
        flat_hash_map<uint64_t, uint64_t> colorsCount;
        uint64_t readID = 0, groupID = 1;
        string seqName, groupName;
        string line;
        priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
        flat_hash_map<string, uint64_t> groupCounter;

        int total_sigs_number = 0;
        int kframe_kSize = selective_kSize;
        if (selective_kSize > 31) kframe_kSize = 31;
        frame = new kDataFramePHMAP(kframe_kSize, mumur_hasher);



        auto* colors = new deduplicatedColumn< StringColorColumn>();
        frame->addColumn("color", (Column*)colors);

        for (const auto& dirEntry : glob2(sigs_dir + "/*")) {
            string file_name = (string)dirEntry;
            size_t lastindex = file_name.find_last_of(".");
            string sig_prefix = file_name.substr(0, lastindex);
            std::string sig_basename = sig_prefix.substr(sig_prefix.find_last_of("/\\") + 1);



            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if (idx != std::string::npos) extension = file_name.substr(idx + 1);
            if (extension != "sig" && extension != "gz") continue;

            zstr::ifstream tmp_stream(file_name);
            JSON sig(tmp_stream);
            int number_of_sub_sigs = sig[0]["signatures"].size();
            string general_name = sig[0]["name"].as<std::string>();
            if (general_name == "") {
                general_name = sig_basename;
            }

            // uncomment if if groupname = sig_name
            // total_sigs_number++;

            for (int i = 0; i < number_of_sub_sigs; i++) {
                int current_kSize = sig[0]["signatures"][i]["ksize"].as<int>();
                if (current_kSize != selective_kSize) continue;

                total_sigs_number++;


                // std::string sig_basename = sig_prefix.substr(sig_prefix.find_last_of("/\\") + 1);
                string md5sum = sig[0]["signatures"][i]["md5sum"].as<std::string>();
                string sig_name = md5sum + ":" + general_name;

                // Here we can decide
                seqName = sig_basename;
                groupName = sig_basename;

                namesMap.insert(make_pair(seqName, groupName));
                auto it = groupNameMap.find(groupName);
                groupCounter[groupName]++;
                if (it == groupNameMap.end()) {
                    groupNameMap.insert(make_pair(groupName, groupID));
                    tagsMap.insert(make_pair(to_string(groupID), groupID));
                    vector<uint32_t> tmp;
                    tmp.clear();
                    tmp.push_back(groupID);
                    legend->insert(make_pair(groupID, tmp));
                    colorsCount.insert(make_pair(groupID, 0));
                    groupID++;
                }
            }
        }

        cout << "namesmap construction done..." << endl;


        flat_hash_map<uint64_t, string> inv_groupNameMap;
        for (auto& _ : groupNameMap)
            inv_groupNameMap[_.second] = _.first;

        string kmer;
        int __batch_count = 0;
        readID = 0;
        int processed_sigs_count = 0;

        for (const auto& dirEntry : glob2(sigs_dir + "/*")) {
            string file_name = (string)dirEntry;
            size_t lastindex = file_name.find_last_of(".");
            string sig_prefix = file_name.substr(0, lastindex);

            std::string sig_basename = sig_prefix.substr(sig_prefix.find_last_of("/\\") + 1);


            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if (idx != std::string::npos) extension = file_name.substr(idx + 1);
            if (extension != "sig" && extension != "gz") continue;


            zstr::ifstream sig_stream(file_name);
            JSON sig(sig_stream);
            int number_of_sub_sigs = sig[0]["signatures"].size();
            string general_name = sig[0]["name"].as<std::string>();
            if (general_name == "") {
                std::string sig_basename = sig_prefix.substr(sig_prefix.find_last_of("/\\") + 1);
                general_name = sig_basename;
            }

            // start

            for (int i = 0; i < number_of_sub_sigs; i++) {

                int current_kSize = sig[0]["signatures"][i]["ksize"].as<int>();
                if (current_kSize != selective_kSize) continue;

                cout << "Processing " << ++processed_sigs_count << "/" << total_sigs_number << " | " << general_name << " k:" << selective_kSize << " ... " << endl;

                string md5sum = sig[0]["signatures"][i]["md5sum"].as<std::string>();
                string sig_name = md5sum + ":" + general_name;
                string readName = sig_basename;
                string groupName = sig_basename;

                flat_hash_map<uint64_t, uint64_t> convertMap;
                auto loaded_sig_it = sig[0]["signatures"][i]["mins"].as_array().begin();
                uint64_t readTag = groupNameMap.find(groupName)->second;
                convertMap.clear();
                convertMap.insert(make_pair(0, readTag));
                convertMap.insert(make_pair(readTag, readTag));

                while (loaded_sig_it != sig[0]["signatures"][i]["mins"].as_array().end()) {
                    string groupName = sig_basename;
                    uint64_t hashed_kmer = loaded_sig_it->as<uint64_t>();

                    frame->insert(hashed_kmer);

                    uint64_t kmerOrder = frame->getkmerOrder(hashed_kmer);
                    if (colors->size() < frame->size() + 1)
                        colors->resize(frame->size() + 1);
                    uint64_t currentTag = colors->index[kmerOrder];
                    auto itc = convertMap.find(currentTag);
                    if (itc == convertMap.end()) {
                        vector<uint32_t> colors = legend->find(currentTag)->second;
                        auto tmpiT = find(colors.begin(), colors.end(), readTag);
                        if (tmpiT == colors.end()) {
                            colors.push_back(readTag);
                            sort(colors.begin(), colors.end());
                        }

                        string colorsString = to_string(colors[0]);
                        for (unsigned int k = 1; k < colors.size(); k++) {
                            colorsString += ";" + to_string(colors[k]);
                        }

                        auto itTag = tagsMap.find(colorsString);
                        if (itTag == tagsMap.end()) {
                            uint64_t newColor;
                            if (freeColors.empty()) {
                                newColor = groupID++;
                            }
                            else {
                                newColor = freeColors.top();
                                freeColors.pop();
                            }

                            tagsMap.insert(make_pair(colorsString, newColor));
                            legend->insert(make_pair(newColor, colors));
                            itTag = tagsMap.find(colorsString);
                            colorsCount[newColor] = 0;
                        }
                        uint64_t newColor = itTag->second;

                        convertMap.insert(make_pair(currentTag, newColor));
                        itc = convertMap.find(currentTag);
                    }

                    if (itc->second != currentTag) {

                        colorsCount[currentTag]--;
                        if (colorsCount[currentTag] == 0 && currentTag != 0) {
                            auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                            if (_invGroupNameIT == inv_groupNameMap.end()) {
                                freeColors.push(currentTag);

                                vector<uint32_t> colors = legend->find(currentTag)->second;
                                string colorsString = to_string(colors[0]);
                                for (unsigned int k = 1; k < colors.size(); k++) {
                                    colorsString += ";" + to_string(colors[k]);
                                }
                                tagsMap.erase(colorsString);

                                legend->erase(currentTag);
                                if (convertMap.find(currentTag) != convertMap.end())
                                    convertMap.erase(currentTag);
                            }

                        }
                        colorsCount[itc->second]++;
                    }

                    colors->index[kmerOrder] = itc->second;


                    readID += 1;
                    groupCounter[groupName]--;
                    if (colorsCount[readTag] == 0) {
                        if (groupCounter[groupName] == 0) {
                            freeColors.push(readTag);
                            legend->erase(readTag);
                        }
                    }
                    loaded_sig_it++;
                }

            }
        }

        colors->values = new StringColorColumn(legend, groupCounter.size());
        delete legend;
        for (auto& iit : namesMap) {
            uint32_t sampleID = groupNameMap[iit.second];
            colors->values->namesMap[sampleID] = iit.second;
        }

        cout << "saving to " << dir_prefix << " ..." << endl;
        frame->save(dir_prefix);

        ofstream f_namesmap;
        f_namesmap.open(dir_prefix + ".kSpider_namesMap");
        for (std::pair<const unsigned int, class std::basic_string<char> >& name : colors->values->namesMap)
            f_namesmap << name.first << "," << name.second << endl;
        f_namesmap.close();

    }

}
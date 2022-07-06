#include "kSpider.hpp"
#include <iostream>
#include <cstdint>
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include "algorithms.hpp"
#include <glob.h>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <string>
#include "RSJparser.tcc"
#include <fstream>

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

    void sourmash_sigs_indexing(string sigs_dir, int selective_kSize) {

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
        int detected_kSize = 0;

        int total_sigs_number = 0;
        frame = new kDataFramePHMAP(selective_kSize, mumur_hasher);

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

            total_sigs_number++;

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

        cout << "namesmap construction done..." << endl;


        // ----------------------------------------------------------------


        flat_hash_map<uint64_t, string> inv_groupNameMap;
        for (auto& _ : groupNameMap)
            inv_groupNameMap[_.second] = _.first;


        int currIndex = 0;
        string kmer;
        uint64_t tagBits = 0;
        uint64_t maxTagValue = (1ULL << tagBits) - 1;

        uint64_t lastTag = 0;
        readID = 0;

        int processed_sigs_count = 0;

        // START
        for (const auto& dirEntry : glob2(sigs_dir + "/*")) {
            string file_name = (string)dirEntry;
            size_t lastindex = file_name.find_last_of(".");
            string sig_prefix = file_name.substr(0, lastindex);

            std::string sig_basename = sig_prefix.substr(sig_prefix.find_last_of("/\\") + 1);

            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if (idx != std::string::npos) extension = file_name.substr(idx + 1);
            if (extension != "sig") continue;

            std::ifstream sig_stream(file_name);
            JSON sig(sig_stream);
            int number_of_sub_sigs = sig[0]["signatures"].size();
            string general_name = sig[0]["name"].as<std::string>();
            if (general_name == "") {
                std::string sig_basename = sig_prefix.substr(sig_prefix.find_last_of("/\\") + 1);
                general_name = sig_basename;
            }


            //START
            for (int i = 0; i < number_of_sub_sigs; i++) {
                int current_kSize = sig[0]["signatures"][i]["ksize"].as<int>();
                if (current_kSize != selective_kSize) continue;

                cout << "Processing " << ++processed_sigs_count << "/" << total_sigs_number << " | " << general_name << " k:" << selective_kSize << " ... " << endl;
                string md5sum = sig[0]["signatures"][i]["md5sum"].as<std::string>();
                string sig_name = md5sum + ":" + general_name;

                flat_hash_map<uint64_t, uint64_t> convertMap;

                string readName = sig_name;
                string groupName = general_name;

                uint64_t readTag = groupNameMap.find(groupName)->second;


                convertMap.clear();
                convertMap.insert(make_pair(0, readTag));
                convertMap.insert(make_pair(readTag, readTag));


                auto loaded_sig_it = sig[0]["signatures"][i]["mins"].as_array().begin();
                while (loaded_sig_it != sig[0]["signatures"][i]["mins"].as_array().end()) {
                    uint64_t hashed_kmer = loaded_sig_it->as<uint64_t>();
                    uint64_t currentTag = frame->getCount(hashed_kmer);
                    auto itc = convertMap.find(currentTag);
                    if (itc == convertMap.end()) {
                        vector<uint32_t> colors = legend->find(currentTag)->second;
                        auto tmpiT = find(colors.begin(), colors.end(), readTag);
                        if (tmpiT == colors.end()) {
                            colors.push_back(readTag);
                            sort(colors.begin(), colors.end());
                        }

                        string colorsString = to_string(colors[0]);
                        for (int k = 1; k < colors.size(); k++) {
                            colorsString += ";" + to_string(colors[k]);
                        }

                        auto itTag = tagsMap.find(colorsString);
                        if (itTag == tagsMap.end()) {
                            uint64_t newColor;
                            if (freeColors.size() == 0) {
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

                    frame->setCount(hashed_kmer, itc->second);
                    if (frame->getCount(hashed_kmer) != itc->second) {
                        //frame->setC(kmer,itc->second);
                        cout << "Error Founded " << hashed_kmer << " from sequence " << readName << " expected "
                            << itc->second << " found " << frame->getCount(hashed_kmer) << endl;
                        exit(1);
                    }
                    loaded_sig_it++;
                }
                readID += 1;
                groupCounter[groupName]--;
                if (colorsCount[readTag] == 0) {
                    if (groupCounter[groupName] == 0) {
                        freeColors.push(readTag);
                        legend->erase(readTag);
                    }

                }
                cout << "   saved_kmers(~" << frame->size() << ")." << endl << endl;
            }
            // END

        }


        colorTable* colors = new intVectorsTable();
        for (auto it : *legend) {
            colors->setColor(it.first, it.second);
        }

        colored_kDataFrame* res = new colored_kDataFrame();
        res->setColorTable(colors);
        res->setkDataFrame(frame);
        for (auto iit = namesMap.begin(); iit != namesMap.end(); iit++) {
            uint32_t sampleID = groupNameMap[iit->second];
            res->namesMap[sampleID] = iit->second;
            res->namesMapInv[iit->second] = sampleID;
        }
        cout << "saving to " << dir_prefix << " ..." << endl;
        res->save(dir_prefix);
    }

}
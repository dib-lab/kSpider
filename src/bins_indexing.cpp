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
#include <fstream>
#include "parallel_hashmap/phmap_dump.h"



using BINS_MAP = phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>,
    phmap::priv::hash_default_hash<std::string>,
    phmap::priv::hash_default_eq<std::string>,
    std::allocator<std::pair<const std::string, phmap::flat_hash_set<uint64_t>>>,
    12,
    std::mutex
>;


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

    void bins_indexing_in_memory(string bins_dir, int selective_kSize) {

        BINS_MAP bin_to_hashes;

        kDataFrame* frame;
        std::string dir_prefix = bins_dir.substr(bins_dir.find_last_of("/\\") + 1);

        flat_hash_map<string, string> namesMap;
        string names_fileName = bins_dir;

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

        int total_bins_number = 0;
        frame = new kDataFramePHMAP(selective_kSize, mumur_hasher);

        flat_hash_map<string, string> basename_to_path;

        for (const auto& dirEntry : glob2(bins_dir + "/*")) {
            string file_name = (string)dirEntry;
            size_t lastindex = file_name.find_last_of(".");
            string bin_prefix = file_name.substr(0, lastindex);
            std::string bin_basename = bin_prefix.substr(bin_prefix.find_last_of("/\\") + 1);


            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if (idx != std::string::npos) extension = file_name.substr(idx + 1);
            if (extension != "bin") {
                cerr << "skipping " << file_name << " does not have extension .bin" << endl;
                continue;
            }

            // Loading bin, insert into namesmap if it contains hashes.
            phmap::flat_hash_set<uint64_t> table_in;
            phmap::BinaryInputArchive ar_in(file_name.c_str());
            table_in.phmap_load(ar_in);
            if (!table_in.size()) continue;

            bin_to_hashes.insert(pair(bin_basename, table_in));
            basename_to_path.insert(pair(bin_basename, file_name));

            total_bins_number++;

            seqName = bin_basename;
            groupName = bin_basename;

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

        int processed_bins_count = 0;

        // START
        for (const auto& [bin_basename, bin_hashes] : bin_to_hashes) {
            //START

            cout << "Processing " << ++processed_bins_count << "/" << total_bins_number << " | " << bin_basename << " ... " << endl;

            flat_hash_map<uint64_t, uint64_t> convertMap;

            string readName = bin_basename;
            string groupName = bin_basename;

            uint64_t readTag = groupNameMap.find(groupName)->second;


            convertMap.clear();
            convertMap.insert(make_pair(0, readTag));
            convertMap.insert(make_pair(readTag, readTag));

            for (const uint64_t& hashed_kmer : bin_hashes) {
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
                            for (int k = 1; k < colors.size(); k++) {
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
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
#include "cpp-json/json.h"
#include "zstr.hpp"
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"

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
        while (sigs_dir.size() > 0 && sigs_dir[sigs_dir.size() - 1] == '/') sigs_dir.erase(sigs_dir.size() - 1, 1);

        std::string dir_prefix = sigs_dir.substr(sigs_dir.find_last_of("/\\") + 1);

        while (dir_prefix.size() > 0 && dir_prefix[dir_prefix.size() - 1] == '/') {
            dir_prefix.erase(dir_prefix.size() - 1, 1);
        }

        cout << "dir_prefix: " << dir_prefix << endl;

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

        flat_hash_map<string, uint32_t> groupName_to_kmerCount;


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
            if (extension != "sig"){
                cerr << "Skipping file: " << file_name << ", must end with .sig even if compressed." << endl;
            };

            zstr::ifstream sig_stream(file_name);
            json::value json = json::parse(sig_stream);


            // select sig based on kSize
            int selected_sig_id = -1;
            for (int sig_no = 0; sig_no < json.size(); sig_no++) {
                json::array& sig_array = as_array(json[sig_no]["signatures"]);
                for (auto it = sig_array.begin(); it != sig_array.end(); ++it) {
                    const json::value& v = *it;
                    if (v["ksize"] == selective_kSize) {
                        selected_sig_id = sig_no;
                        break;
                    }
                }
            }
            if (selected_sig_id == -1) {
                throw std::runtime_error("No signature found for kSize: " + std::to_string(selective_kSize));
            }


            auto sourmash_sig = json[selected_sig_id]["signatures"];
            const json::array& sig_array = as_array(sourmash_sig);


            //START
            for (auto it = sig_array.begin(); it != sig_array.end(); ++it) {

                const json::value& v = *it;
                if (v["ksize"] != selective_kSize) {
                    continue;
                }

                cout << "Processing " << ++processed_sigs_count << "/" << total_sigs_number << " | " << sig_basename << " k:" << selective_kSize << " ... " << endl;


                flat_hash_map<uint64_t, uint64_t> convertMap;

                string readName = sig_basename;
                string groupName = sig_basename;

                uint64_t readTag = groupNameMap.find(groupName)->second;


                convertMap.clear();
                convertMap.insert(make_pair(0, readTag));
                convertMap.insert(make_pair(readTag, readTag));

                const json::array& mins = as_array(v["mins"]);
                auto loaded_sig_it = mins.begin();
                groupName_to_kmerCount[groupName] = mins.size();


                while (loaded_sig_it != mins.end()) {
                    uint64_t hashed_kmer = json::to_number<uint64_t>(*loaded_sig_it);
                    uint64_t currentTag = frame->getCount(hashed_kmer);
                    auto itc = convertMap.find(currentTag);
                    if (itc == convertMap.end()) {
                        vector<uint32_t> colors = legend->find(currentTag)->second;
                        auto tmpiT = find(colors.begin(), colors.end(), readTag);
                        if (tmpiT == colors.end()) {
                            colors.push_back(readTag);
                            sort(colors.begin(), colors.end());
                        }


                        // [OPTIMIZE] [TODO] Optimize colors concatenation
                        // string colorsString = to_string(colors[0]);
                        // for (int k = 1; k < colors.size(); k++) {
                        //     colorsString += ";" + to_string(colors[k]);
                        // }

                        std::stringstream ss;
                        if (!colors.empty()) {
                            ss << colors[0];
                            for (size_t k = 1; k < colors.size(); ++k) {
                                ss << ";" << colors[k];
                            }
                        }
                        std::string colorsString = ss.str();
                        // END [OPTIMIZE] [TODO] Optimize colors concatenation


                        auto itTag = tagsMap.find(colorsString);

                        // START [OPTIMIZE] [TODO] Optimize
                        // if (itTag == tagsMap.end()) {
                        //     uint64_t newColor;
                        //     if (freeColors.size() == 0) {
                        //         newColor = groupID++;
                        //     }
                        //     else {
                        //         newColor = freeColors.top();
                        //         freeColors.pop();
                        //     }
                        //     tagsMap.insert(make_pair(colorsString, newColor));
                        //     legend->insert(make_pair(newColor, colors));
                        //     itTag = tagsMap.find(colorsString);
                        //     colorsCount[newColor] = 0;
                        // }

                        if (itTag == tagsMap.end()) {
                            uint64_t newColor;
                            if (freeColors.empty()) newColor = groupID++;
                            else {
                                newColor = freeColors.top();
                                freeColors.pop();
                            }
                            auto inserted = tagsMap.emplace(colorsString, newColor);
                            legend->emplace(newColor, colors);
                            itTag = inserted.first;
                            colorsCount[newColor] = 0;
                        }
                        // END [OPTIMIZE] [TODO] Optimize
                        uint64_t newColor = itTag->second;

                        convertMap.emplace(currentTag, newColor);
                        itc = convertMap.find(currentTag);
                    }
                    // START [OPTIMIZE] [TODO] Optimize
                    // if (itc->second != currentTag) {

                    //     colorsCount[currentTag]--;
                    //     if (colorsCount[currentTag] == 0 && currentTag != 0) {

                    //         auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                    //         if (_invGroupNameIT == inv_groupNameMap.end()) {
                    //             freeColors.push(currentTag);
                    //             vector<uint32_t> colors = legend->find(currentTag)->second;
                    //             string colorsString = to_string(colors[0]);
                    //             for (unsigned int k = 1; k < colors.size(); k++) {
                    //                 colorsString += ";" + to_string(colors[k]);
                    //             }
                    //             tagsMap.erase(colorsString);
                    //             legend->erase(currentTag);
                    //             if (convertMap.find(currentTag) != convertMap.end())
                    //                 convertMap.erase(currentTag);
                    //         }

                    //     }
                    //     colorsCount[itc->second]++;
                    // }

                    if (itc->second != currentTag) {
                        --colorsCount[currentTag];
                        if (colorsCount[currentTag] == 0 && currentTag != 0) {
                            auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                            if (_invGroupNameIT == inv_groupNameMap.end()) {
                                freeColors.push(currentTag);
                                auto legendIt = legend->find(currentTag);
                                std::vector<uint32_t> colors = legendIt->second;
                                std::stringstream ss;
                                ss << colors[0];
                                for (unsigned int k = 1; k < colors.size(); k++) {
                                    ss << ";" << colors[k];
                                }
                                std::string colorsString = ss.str();
                                tagsMap.erase(colorsString);
                                legend->erase(legendIt);
                                auto convertMapIt = convertMap.find(currentTag);
                                if (convertMapIt != convertMap.end()) {
                                    convertMap.erase(convertMapIt);
                                }
                            }
                        }
                        ++colorsCount[itc->second];
                    }

                    // END [OPTIMIZE] [TODO] Optimize


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
                cout << "   saved_kmers(~" << frame->size() << ")." << endl;
                cout << "   colors(~" << legend->size() << ")." << endl << endl;

                break;
            }
            // END

        }

        cerr << "Indexing done..." << endl;
        cerr << "Total number of colors: " << legend->size() << endl;


        string output_prefix = dir_prefix;

        // Dump kmer count
        flat_hash_map<uint32_t, uint32_t> groupID_to_kmerCount;
        for (auto& [groupName, kmerCount] : groupName_to_kmerCount) {
            groupID_to_kmerCount[groupNameMap[groupName]] = kmerCount;
        }

        phmap::BinaryOutputArchive ar_out(string(output_prefix + "_groupID_to_kmerCount.bin").c_str());
        groupID_to_kmerCount.phmap_dump(ar_out);

        // write legend in csv
        ofstream legend_file;
        legend_file.open(output_prefix + "_legend.csv");
        legend_file << "color,groupIDs" << endl;
        cout << "Legend size: " << legend->size() << endl;
        for (auto it : *legend) {
            legend_file << it.first << ",";
            for (auto it2 : it.second) {
                legend_file << it2 << ";";
            }
            legend_file << endl;
        }
        legend_file.close();

        // Write color sizes
        ofstream color_sizes_file;
        color_sizes_file.open(output_prefix + "_color_sizes.csv");
        color_sizes_file << "color,size" << endl;
        for (auto it : *legend) {
            color_sizes_file << it.first << "," << it.second.size() << endl;
        }
        color_sizes_file.close();

        // write color counts
        ofstream color_counts_file;
        color_counts_file.open(output_prefix + "_color_counts.csv");
        color_counts_file << "color,count" << endl;
        for (auto it : colorsCount) {
            color_counts_file << it.first << "," << it.second << endl;
        }
        color_counts_file.close();
        


       // Dump color->sources
        double average_color_size = 0;
        double color_size_standard_deviation = 0;
        auto color_to_sources = new phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint32_t>>();
        for (auto it : *legend) {
            if (colorsCount[it.first] == 0) continue;
            color_size_standard_deviation += pow(it.second.size() - average_color_size, 2);
            average_color_size += it.second.size();
            phmap::flat_hash_set<uint32_t> tmp(std::make_move_iterator(it.second.begin()), std::make_move_iterator(it.second.end()));
            color_to_sources->operator[](it.first) = tmp;
        }
        average_color_size /= legend->size();
        color_size_standard_deviation /= legend->size();
        color_size_standard_deviation = sqrt(color_size_standard_deviation);
        cerr << "Total selected colors: " << color_to_sources->size() << endl;
        cerr << "Average color size: " << average_color_size << endl;
        cerr << "Color size standard deviation: " << color_size_standard_deviation << endl;


        phmap::BinaryOutputArchive ar_out_1(string(output_prefix + "_color_to_sources.bin").c_str());
        ar_out_1.saveBinary(color_to_sources->size());
        for (auto& [k, v] : *color_to_sources)
        {
            ar_out_1.saveBinary(k);
            ar_out_1.saveBinary(v);
        }

        // Dump colors count
        phmap::BinaryOutputArchive ar_out_3(string(output_prefix + "_color_count.bin").c_str());
        colorsCount.phmap_dump(ar_out_3);

        // export namesMap
        ofstream namesMapOut(output_prefix + ".namesMap");
        namesMapOut << namesMap.size() << endl;
        for (auto it : namesMap)
        {
            namesMapOut << groupNameMap[it.second] << " " << it.second << endl;
        }
        namesMapOut.close();

        // Write extra info
        ofstream file(output_prefix + ".extra");
        file << selective_kSize << endl;
        file << frame->KD->hash_mode << endl;
        file << frame->KD->slicing_mode << endl;
        file << frame->KD->params_to_string() << endl;
        file << "features:" << frame->getkSize() << endl;
        file.close();

        // ------- Pause serializing index for now.
        /*

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
        */
        // ------ END Pause serializing index for now.

    }

}
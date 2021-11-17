#include "kSpider.hpp"
#include <iostream>
#include <cstdint>
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include "algorithms.hpp"
#include <filesystem>
using recursive_directory_iterator = std::filesystem::recursive_directory_iterator;


namespace kSpider {

    void index_kmers(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix) {
        auto* KF = new kDataFrameMQF(KMERS, integer_hasher, { {"kSize", kSize} });
        auto* ckf = kProcessor::index(KF, fasta_file, chunk_size, names_file);
        ckf->save(index_prefix);
    }

    void index_kmers_nonCanonical(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix) {
        auto* KF = new kDataFrameMQF(KMERS, nonCanonicalInteger_Hasher, { {"kSize", kSize} });
        auto* ckf = kProcessor::index(KF, fasta_file, chunk_size, names_file);
        ckf->save(index_prefix);
    }

    void index_skipmers(int m, int n, int k, string fasta_file, string names_file, int chunk_size, string index_prefix) {
        auto* KF = new kDataFramePHMAP(SKIPMERS, integer_hasher, { {"m", m}, {"n", n}, {"k", k} });
        auto* ckf = kProcessor::index(KF, fasta_file, chunk_size, names_file);
        ckf->save(index_prefix);
    }

    void index_protein(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix) {
        auto* KF = new kDataFramePHMAP(PROTEIN, protein_hasher, { {"kSize", kSize} });
        auto* ckf = kProcessor::index(KF, fasta_file, chunk_size, names_file);
        ckf->save(index_prefix);
    }

    void index_dayhoff(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix) {
        auto* KF = new kDataFramePHMAP(PROTEIN, proteinDayhoff_hasher, { {"kSize", kSize} });
        auto* ckf = kProcessor::index(KF, fasta_file, chunk_size, names_file);
        ckf->save(index_prefix);
    }


    void index_datasets(string kfs_dir) {

        kDataFrame * frame;
        
        
        std::filesystem::path dirname(kfs_dir);
        string dir_prefix = dirname.parent_path().filename();

        for (const auto& dirEntry : recursive_directory_iterator(kfs_dir)) {
            string file_name = (string)dirEntry.path();
            size_t lastindex = file_name.find_last_of(".");
            string kf_prefix = file_name.substr(0, lastindex);
            auto * _kf = kDataFrame::load(kf_prefix);
            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if(idx != std::string::npos) extension = file_name.substr(idx+1);

            if(extension == "mqf") frame = new kDataFrameMQF(_kf->getkSize());
            else if (extension == "phmap") frame = new kDataFramePHMAP(_kf->getkSize());

        }
        

        flat_hash_map<string, string> namesMap;
        flat_hash_map<string, uint64_t> tagsMap;
        flat_hash_map<string, uint64_t> groupNameMap;
        auto* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
        flat_hash_map<uint64_t, uint64_t> colorsCount;
        uint64_t readID = 0, groupID = 1;
        string seqName, groupName;
        string line;
        priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
        flat_hash_map<string, uint64_t> groupCounter;

        flat_hash_map<uint64_t, string> inv_groupNameMap;
        for (auto& _ : groupNameMap)
            inv_groupNameMap[_.second] = _.first;


        vector<kDataFrameMQF*> frames;
        int currIndex = 0;
        string kmer;
        uint64_t tagBits = 0;
        uint64_t maxTagValue = (1ULL << tagBits) - 1;
        //  kDataFrame *frame;

        uint64_t lastTag = 0;
        readID = 0;
        int __batch_count = 0;

        for (const auto& dirEntry : recursive_directory_iterator(kfs_dir)) {

            string file_name = (string)dirEntry.path();
            size_t lastindex = file_name.find_last_of(".");
            string kf_prefix = file_name.substr(0, lastindex);
            auto * loaded_kf = kDataFrame::load(kf_prefix);
            

            flat_hash_map<uint64_t, uint64_t> convertMap;
            string readName = kf_prefix;
            auto loaded_kf_it = loaded_kf->begin();

                auto it = namesMap.find(readName);
                if (it == namesMap.end()) {
                    continue;
                    // cout << "read " << readName << "dont have group. Please check the group names file." << endl;
                }
                string groupName = it->second;

                uint64_t readTag = groupNameMap.find(groupName)->second;


                convertMap.clear();
                convertMap.insert(make_pair(0, readTag));
                convertMap.insert(make_pair(readTag, readTag));

               while (loaded_kf_it != loaded_kf->end()) {
                    uint64_t hashed_kmer = loaded_kf_it.getHashedKmer();
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
                            // if(groupID>=maxTagValue){
                            //   cerr<<"Tag size is not enough. ids reached "<<groupID<<endl;
                            //   return -1;
                            // }
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
                    it++;
                }
                readID += 1;
                groupCounter[groupName]--;
                if (colorsCount[readTag] == 0) {
                    if (groupCounter[groupName] == 0) {
                        freeColors.push(readTag);
                        legend->erase(readTag);
                    }
                    
                }
                
            
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

        res->save(dir_prefix);
    }




}
#include "kSpider.hpp"
#include <iostream>
#include <cstdint>
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include "algorithms.hpp"
#include <glob.h> // glob(), globfree()
#include <vector>
#include <stdexcept>
#include <sstream>
#include <string.h>
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"

#define LOWEST_PERCENTILE 5


// thanks to https://stackoverflow.com/a/8615450/3371177
std::vector<std::string> glob(const std::string& pattern) {
    using namespace std;

    // glob struct resides on the stack
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // do the glob operation
    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if(return_value != 0) {
        globfree(&glob_result);
        stringstream ss;
        ss << "glob() failed with return_value " << return_value << endl;
        throw std::runtime_error(ss.str());
    }

    // collect all the filenames into a std::list<std::string>
    vector<string> filenames;
    for(size_t i = 0; i < glob_result.gl_pathc; ++i) {
        filenames.push_back(string(glob_result.gl_pathv[i]));
    }

    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}


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
        std::string dir_prefix = kfs_dir.substr(kfs_dir.find_last_of("/\\") + 1);

        flat_hash_map<string, string> namesMap;
        string names_fileName = kfs_dir;

        flat_hash_map<string, uint64_t> tagsMap;
        flat_hash_map<string, uint64_t> groupNameMap;
        flat_hash_map<uint64_t, std::vector<uint32_t>> *legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
        flat_hash_map<uint64_t, uint64_t> colorsCount;
        uint64_t readID = 0, groupID = 1;
        string seqName, groupName;
        string line;
        priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
        flat_hash_map<string, uint64_t> groupCounter;

        int total_kfs_number = 0;
        int detected_kSize;

        // get kSize and type
         for (const auto& dirEntry : glob(kfs_dir + "/*")) {
            string file_name = (string)dirEntry;
            size_t lastindex = file_name.find_last_of(".");
            string kf_prefix = file_name.substr(0, lastindex);
            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if(idx != std::string::npos) extension = file_name.substr(idx+1);
            
            hashingModes _hm;
            if(extension == "mqf" || extension == "phmap") {
                auto * _kf = kDataFrame::load(kf_prefix);
                _hm = _kf->getkmerDecoder()->hash_mode;
                detected_kSize = _kf->getkSize();
                cout << "Detected kSize: " << detected_kSize << endl;
            }else{
                continue;
            }
            // if(extension == "mqf") {frame = new kDataFrameMQF(detected_kSize, 30, mumur_hasher); break;} // temp. switch off
            if(extension == "mqf") {frame = new kDataFramePHMAP(detected_kSize, _hm); break;}
            else if(extension == "phmap") {frame = new kDataFramePHMAP(detected_kSize, _hm); break;}
            else {continue;}
        }

        for (const auto& dirEntry : glob(kfs_dir + "/*")) {
            string file_name = (string)dirEntry;
            size_t lastindex = file_name.find_last_of(".");
            string kf_prefix = file_name.substr(0, lastindex);

            
            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if(idx != std::string::npos) extension = file_name.substr(idx+1);
            if(extension != "mqf" and extension != "phmap") continue;

            total_kfs_number++;

            std::string kf_basename = kf_prefix.substr(kf_prefix.find_last_of("/\\") + 1);

            seqName = kf_basename;
            groupName = kf_basename;
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
        //  kDataFrame *frame;

        uint64_t lastTag = 0;
        readID = 0;
        int __batch_count = 0;

        int processed_kfs_count = 0;
        flat_hash_map<string, uint32_t> groupName_to_kmerCount;

        for (const auto& dirEntry : glob(kfs_dir + "/*")) {
            string file_name = (string)dirEntry;
            size_t lastindex = file_name.find_last_of(".");
            string kf_prefix = file_name.substr(0, lastindex);

            std::string kf_basename = kf_prefix.substr(kf_prefix.find_last_of("/\\") + 1);
            
            
            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if(idx != std::string::npos) extension = file_name.substr(idx+1);
            if(extension != "mqf" and extension != "phmap") continue;
            
            // uint64_t removed_kmers_from_percentile = 0;
            auto * loaded_kf = kDataFrame::load(kf_prefix);
            cout << "Processing " << ++processed_kfs_count << "/" << total_kfs_number << " | " << kf_basename << " k:" << loaded_kf->ksize() << " ... " << endl;

            // Calculating percentile
            // auto perc_kf_it = loaded_kf->begin();
            // auto* kmerCounts = new(vector<uint64_t>);
            // uint64_t real_no_of_kmers = 0;
            // while (perc_kf_it != loaded_kf->end()) {
            //     kmerCounts->push_back(perc_kf_it.getCount());
            //     real_no_of_kmers++;
            //     perc_kf_it++;
            // }
            // sort(kmerCounts->begin(), kmerCounts->end());
            // uint64_t _idx = (uint64_t)ceil((real_no_of_kmers * LOWEST_PERCENTILE / 100));
            // uint64_t count_percentile_cutoff = kmerCounts->at(_idx);
            // cout << "   calculated percentile cutoff kmercount=" << count_percentile_cutoff << endl;
            // delete kmerCounts;

            flat_hash_map<uint64_t, uint64_t> convertMap;
            string readName = kf_basename;
            auto loaded_kf_it = loaded_kf->begin();

                string groupName = kf_basename;
                groupName_to_kmerCount[groupName] = loaded_kf->size();

                uint64_t readTag = groupNameMap.find(groupName)->second;


                convertMap.clear();
                convertMap.insert(make_pair(0, readTag));
                convertMap.insert(make_pair(readTag, readTag));

               while (loaded_kf_it != loaded_kf->end()) {
                    
                    // Don't consider the kmer if it exist in the lower percentile
                    // if(loaded_kf_it.getCount() < count_percentile_cutoff){
                    //     removed_kmers_from_percentile++;
                    //     loaded_kf_it++;
                    //     continue; 
                    // }

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
                    loaded_kf_it++;
                }
                // cout << "   Removed kmers from percentile=(" << removed_kmers_from_percentile <<") out of ("<< real_no_of_kmers <<")" << endl;
                readID += 1;
                groupCounter[groupName]--;
                if (colorsCount[readTag] == 0) {
                    if (groupCounter[groupName] == 0) {
                        freeColors.push(readTag);
                        legend->erase(readTag);
                    }
                    
                }
                cout << "   saved_kmers(~" << frame->size() << ")." << endl << endl;
                delete loaded_kf;
        }

        string output_prefix = dir_prefix;
        
        // Dump kmer count
        flat_hash_map<uint32_t, uint32_t> groupID_to_kmerCount;
        for(auto & [groupName, kmerCount] : groupName_to_kmerCount){
            groupID_to_kmerCount[groupNameMap[groupName]] = kmerCount;
        }

        phmap::BinaryOutputArchive ar_out(string(output_prefix + "_groupID_to_kmerCount.bin").c_str());
        groupID_to_kmerCount.phmap_dump(ar_out);


        // Dump color->sources

        auto color_to_sources = new phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint32_t>>();
        for (auto it : *legend) {
            phmap::flat_hash_set<uint32_t> tmp(std::make_move_iterator(it.second.begin()), std::make_move_iterator(it.second.end()));
            color_to_sources->operator[](it.first) = tmp;
        }

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


        colorTable* colors = new intVectorsTable();
        for (auto it : *legend) {
            colors->setColor(it.first, it.second);
        }

        // export namesMap
        ofstream namesMapOut(output_prefix + ".namesMap");
        namesMapOut<<namesMap.size()<<endl;
        for(auto it:namesMap)
        {
            namesMapOut<<groupNameMap[it.second]<<" "<<it.second<<endl;
        }
        namesMapOut.close();

        // Write extra info
        ofstream file(output_prefix + ".extra");
        file << detected_kSize << endl;
        file << frame->KD->hash_mode << endl;
        file << frame->KD->slicing_mode << endl;
        file << frame->KD->params_to_string() << endl;
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
        cout << "saving to "<< dir_prefix << " ..." << endl;
        res->save(dir_prefix);
        */
        // ------ END Pause serializing index for now.
    }

}
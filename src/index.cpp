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

        // get kSize and type
         for (const auto& dirEntry : glob(kfs_dir + "/*")) {
            string file_name = (string)dirEntry;
            size_t lastindex = file_name.find_last_of(".");
            string kf_prefix = file_name.substr(0, lastindex);
            std::string::size_type idx;
            idx = file_name.rfind('.');
            std::string extension = "";
            if(idx != std::string::npos) extension = file_name.substr(idx+1);
            if(extension == "mqf") {frame = new kDataFrameMQF(kDataFrame::load(kf_prefix)->getkSize()); break;}
            else if(extension == "mqf") {frame = new kDataFramePHMAP(kDataFrame::load(kf_prefix)->getkSize()); break;}
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



        // ----------------------------------------------------------------


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

        int processed_kfs_count = 0;

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
            
            auto * loaded_kf = kDataFrame::load(kf_prefix);
            cout << "Processing " << ++processed_kfs_count << "/" << total_kfs_number << " | " << kf_basename << " k:" << loaded_kf->ksize() << " ..." << endl;

            flat_hash_map<uint64_t, uint64_t> convertMap;
            string readName = kf_basename;
            auto loaded_kf_it = loaded_kf->begin();

                string groupName = kf_basename;

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
                    loaded_kf_it++;
                }
                readID += 1;
                groupCounter[groupName]--;
                if (colorsCount[readTag] == 0) {
                    if (groupCounter[groupName] == 0) {
                        freeColors.push(readTag);
                        legend->erase(readTag);
                    }
                    
                }
                delete loaded_kf;
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
        cout << "saving to "<< dir_prefix << " ..." << endl;
        res->save(dir_prefix);
    }

}
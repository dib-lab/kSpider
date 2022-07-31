#include <iostream>
#include <cstdint>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include <ctime>
#include <omp.h>
#include "cpp-json/json.h"
#include "zstr.hpp"
#include <glob.h>
#include <string>
#include <stdexcept>
#include "parallel_hashmap/phmap_dump.h"
#include <cstdlib>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

typedef std::chrono::high_resolution_clock Time;

inline bool file_exists(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}


std::vector<std::string> glob2(const std::string& pattern) {
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

    vector<string> filenames;
    for (size_t i = 0; i < glob_result.gl_pathc; ++i)
        filenames.push_back(string(glob_result.gl_pathv[i]));

    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}


int main(int argc, char** argv) {

    if (argc != 5) {
        cout << "run: ./sigs_to_bins <sigs_directory> <kSize> <output_dir> <threads>" << endl;
        exit(1);
    }
    string sigs_dir = argv[1];
    int kSize = stoi(argv[2]);
    string output_dir = argv[3];
    int user_threads = stoi(argv[4]);

    string cmd = "mkdir -p " + output_dir;

    const int dir_err = system(cmd.c_str());
    if (-1 == dir_err)
    {
        printf("Error creating directory!n");
        exit(1);
    }


    // 1. Scan all sigs in a directory
    vector<string> sigs_paths;
    vector<string> sig_names;

    int skipped_files = 0;

    int total_sigs_number = 0;
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

        if (file_exists(output_dir + "/" + sig_basename + ".bin")) { skipped_files++; continue; }

        sig_names.push_back(sig_basename);
        sigs_paths.push_back(file_name);

        total_sigs_number++;
    }

    cout << "Skipped " << skipped_files << " files as they already converted to bins." << endl;

    int sigs_count = sigs_paths.size();
    auto begin_time = Time::now();

#pragma omp parallel num_threads(user_threads)
    {
#pragma omp for
        for (int j = 0; j < sigs_paths.size(); j++) {
            string& sig_path = sigs_paths[j];
            string& sig_name = sig_names[j];
            phmap::flat_hash_set<uint64_t> tmp_hashes;

            zstr::ifstream sig_stream(sig_path);
            json::value json = json::parse(sig_stream);
            auto sourmash_sig = json[0]["signatures"];
            const json::array& sig_array = as_array(sourmash_sig);
            for (auto it = sig_array.begin(); it != sig_array.end(); ++it) {
                const json::value& v = *it;
                if (v["ksize"] == kSize) {
                    const json::array& mins = as_array(v["mins"]);
                    const json::array& abundances = as_array(v["abundances"]);
                    auto mins_it = mins.begin();
                    while (mins_it != mins.end()) {
                        // const auto & abund = json::to_number<int>(*abund_it);
                        tmp_hashes.insert(json::to_number<uint64_t>(*mins_it));
                        mins_it++;
                    }
                }
                break;
            }

            string out_path = output_dir + "/" + sig_name + ".bin";
            phmap::BinaryOutputArchive ar_out(out_path.c_str());
            tmp_hashes.phmap_dump(ar_out);
        }
    }

    cout << endl;
    cout << "Process completed in " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;

}
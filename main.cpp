#include "virtualQs.hpp"
#include <iostream>
#include "stdint.h"
#include "argh.h"
#include "combinations.hpp"
#include <chrono>
#include "pivote.hpp"

typedef std::chrono::high_resolution_clock Time;

int main(int argc, char **argv) {


    // TODO: Finish the aguments parsing later.

    // Arguments parsing ----------------------------

    int minQ, maxQ, stepQ;
    string index_prefix;
    string all_Qs = "";
    set<int> Qs_set;

    argh::parser cmdl(argv);
    bool entered_minQ = (!!(cmdl({"-m", "--min-q"}) >> minQ));
    bool entered_maxQ = (!!(cmdl({"-M", "--max-q"}) >> maxQ));
    bool entered_stepQ = (!!(cmdl({"-s", "--step-q"}) >> stepQ));
    bool entered_index_prefix = (!!(cmdl({"-i", "--idx"}) >> index_prefix));
    bool entered_manualQs = (!!(cmdl({"-l", "--qs"}) >> all_Qs));

    if (!entered_index_prefix) {
        cerr << "must provide an index prefix" << endl;
        exit(1);
    }

    if (entered_manualQs) {
        std::stringstream ss(all_Qs);
        for (int i; ss >> i;) {
            Qs_set.insert(i);
            if (ss.peek() == ',')
                ss.ignore();
        }
    } else if (entered_minQ && entered_maxQ && entered_stepQ) {
        for (int Q = minQ; Q <= maxQ; Q += stepQ)
            Qs_set.insert(Q);
    }

    // End Arguments parsing ----------------------------

    if(cmdl[1] == "pivote"){
        cerr << "pivoting..." << endl;
        pivote(index_prefix, Qs_set);
        exit(0);
    }

    // Instantiating VirtualQs class
    auto t1 = Time::now();
    virtualQs VQ = virtualQs(index_prefix, Qs_set);

    auto t2 = Time::now();
    cerr << "[SUCCESS] Done initializing the virtualQs in: "
         << std::chrono::duration_cast<chrono::seconds>(t2 - t1).count()
         << " secs." << endl;


    cerr << "[INFO] scanning virtualQs ..." << endl;
    t1 = Time::now();

    for (int const &Q : Qs_set) {
        uint64_t mask = VQ.masks[Q];
        VQ.curr_Q = Q;
        auto it = VQ.KF->begin();
        uint64_t prev_kmer = it.getHashedKmer();
        uint64_t prev_kmer_color = it.getKmerCount();
        uint64_t XOR;
        uint64_t curr_kmer = 0;
        uint64_t curr_kmer_color = 0;

        while (it != VQ.KF->end()) {
            it++;
            curr_kmer = it.getHashedKmer();
            curr_kmer_color = it.getKmerCount();
            XOR = prev_kmer xor curr_kmer;

            bool matched = !(bool) (XOR & mask);

            if (matched) {
                VQ.temp_superColors.push_back(prev_kmer_color);
                VQ.temp_superColors.push_back(curr_kmer_color);
            } else {
                VQ.temp_superColors.push_back(prev_kmer_color);
                uint64_t super_color_id = VQ.create_super_color(VQ.temp_superColors);
                bool super_color_exist = (VQ.superColors.find(super_color_id) != VQ.superColors.end());

                if (super_color_exist) {
                    VQ.superColorsCount[super_color_id]++;
                } else {
                    VQ.superColors[super_color_id] = VQ.temp_superColors;
                    VQ.superColorsCount[super_color_id] = 1;
                }

                VQ.temp_superColors.clear();
                VQ.temp_superColors.push_back(curr_kmer_color);
            }
            prev_kmer = curr_kmer;
            prev_kmer_color = curr_kmer_color;
        }


        for (auto &superColor : VQ.temp_superColors) {
            VQ.temp_superColors.erase(std::remove(VQ.temp_superColors.begin(), VQ.temp_superColors.end(), curr_kmer_color), VQ.temp_superColors.end());
            if (VQ.temp_superColors.empty()) {
                continue;
            }

            uint64_t super_color_id = VQ.create_super_color(VQ.temp_superColors);
            bool super_color_exist = (VQ.superColors.find(super_color_id) != VQ.superColors.end());

            if (super_color_exist) {
                VQ.superColorsCount[super_color_id]++;
            } else {
                VQ.superColors[super_color_id] = VQ.temp_superColors;
                VQ.superColorsCount[super_color_id] = 1;
            }

        }

        cerr << "superColors size: " << VQ.superColors.size() << endl;

//         Clearing
        VQ.pairwise();

        VQ.export_to_tsv();
//        VQ.export_to_sqlite();
        VQ.superColors.clear();
//        VQ.edges.clear();

        for (uint64_t i = 0; i < VQ.no_seqs; i++) {
            VQ.edges2[i].clear();
        }

        VQ.superColorsCount.clear();
        VQ.temp_superColors.clear();
    }


//    VQ.SQL->close();
    VQ.myfile.close();
    return EXIT_SUCCESS;
}



//    t2 = Time::now();
//    cerr << "[SUCCESS] Done scanning virtualQs in: "
//         << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
//         << " secs." << endl;
//
//    cerr << "[INFO] constructing pairwise matrix" << endl;
//    t1 = Time::now();
//    VQ.pairwise();
//    t2 = Time::now();
//    cerr << "[SUCCESS] Done constructing pairwise matrix in: "
//         << std::chrono::duration_cast<chrono::seconds>(t2 - t1).count() << " secs." << endl;
//
//    cerr << "[INFO] Exporting to TSV" << endl;
//    t1 = Time::now();
//    VQ.export_to_tsv();
//    t2 = Time::now();
//    cerr << "[SUCCESS] Done Exporting to TSV in: "
//         << std::chrono::duration_cast<chrono::seconds>(t2 - t1).count() << " secs." << endl;

//    for (auto const &edge : VQ.edges) {
//        cout << "seq(" << edge.first.first << "," << edge.first.second << "): ";
//        for(auto const & Q:edge.second){
//            cout << "Q" << (int) Q.first << ":" << Q.second << " | ";
//        }
//        cout << endl;
//
//    }






//    cerr << "[INFO] Combinations.." << endl;
//    Combo combo = Combo();
//    combo.combinations(VQ.seq_to_kmers_no.size());

//    for (auto const &seq_pair : combo.combs) {
//        if(seq_pair.first < seq_pair.second)
//            continue;
//
//        for(auto const & Q : VQ.mainQs){
//            uint32_t _seq1 = seq_pair.first + 1;
//            uint32_t _seq2 = seq_pair.second + 1;
//            cout << "seq(" << _seq1 << "," << _seq2 << "): ";
//            cout << "Q" << (int)Q << " = ";
//            cout << VQ.edges[{{_seq1, _seq2},Q}];
//            cout << endl;
//        }
//    }

//    for (auto const &edge : VQ.edges) {
//        cout << "seq(" << edge.first.first.first << "," << edge.first.first.second << "): ";
//        cout << "Q" << (int)edge.first.second << " = " << edge.second << endl;
//    }

//    VQ.export_to_sqlite();

//    VQ.fast_export_to_tsv();
//    cerr << "Exporting" << endl;
//    VQ.fast_export_to_tsv();

//    for (auto const &edge : VQ.edges) {
//        cout << "seq(" << edge.first.first.first << "," << edge.first.first.second << "): ";
//        cout << "Q" << (int)edge.first.second << " = " << edge.second << endl;
//    }



//     // Exporting superColors

//    // for (auto &superColor : VQ.superColors) {
//     int Q = superColor.first;
//     cout << "Q" << Q << endl;
//     for (auto const &color : superColor.second) {
//         uint64_t supercolor_id = color.first;
//         cout << supercolor_id << ": ";
//         for (auto const &c : color.second) {
//             cout << c << ", ";
//         }
//         cout << endl;
//     }
//     cout << endl;
// }
// cout << "---------------\nCount: \n";

//    // for (auto &superColor : VQ.superColorsCount) {
//     int Q = superColor.first;
//     cout << "Q" << Q << endl;
//     for (auto const &color : superColor.second) {
//         uint64_t supercolor_id = color.first;
//         cout << supercolor_id << ": " << color.second;
//         cout << endl;
//     }
//     cout << endl;
// }
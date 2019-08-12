#include "virtualQs.hpp"
#include <iostream>
#include "stdint.h"
#include "argh.h"

int main(int argc, char **argv) {

    // TODO: Finish the aguments parsing later.
    // TODO: Change the Qs sys args input for excluding written Qs.

    int minQ, maxQ, stepQ;
    string index_prefix;

    argh::parser cmdl(argv);

    cmdl({"-m", "--min-q"}) >> minQ;
    cmdl({"-M", "--max-q"}) >> maxQ;
    cmdl({"-s", "--step-q"}) >> stepQ;
    cmdl({"-i", "--idx"}) >> index_prefix;

    // Arguments parsing ----------------------------

    virtualQs VQ = virtualQs(index_prefix, minQ, maxQ, stepQ);

    auto it = VQ.KF->begin();

    uint64_t prev_kmer = it.getHashedKmer();
    uint64_t prev_kmer_color = it.getKmerCount();
    uint64_t XOR;
    uint64_t curr_kmer;
    uint64_t curr_kmer_color;

    bool matched;

    while (it != VQ.KF->end()) {
        it++;
        curr_kmer = it.getHashedKmer();
        curr_kmer_color = it.getKmerCount();
        XOR = prev_kmer xor curr_kmer;

        for (auto const &mask : VQ.masks) {
            int Q = mask.first;
            matched = !(bool) (XOR & mask.second);

            if (matched) {
                VQ.temp_superColors[Q].insert(prev_kmer_color);
                VQ.temp_superColors[Q].insert(curr_kmer_color);
            } else {
                VQ.temp_superColors[Q].insert(prev_kmer_color);
                uint64_t super_color_id = VQ.create_super_color(VQ.temp_superColors[Q]);
                bool super_color_exist = (VQ.superColors[Q].find(super_color_id) != VQ.superColors[Q].end());

                if (super_color_exist) {
                    VQ.superColorsCount[Q][super_color_id]++;
                } else {
                    VQ.superColors[Q][super_color_id] = VQ.temp_superColors[Q];
                    VQ.superColorsCount[Q][super_color_id] = 1;
                }

                VQ.temp_superColors[Q].clear();
                VQ.temp_superColors[Q].insert(curr_kmer_color);
            }

        }


        prev_kmer = curr_kmer;
        prev_kmer_color = curr_kmer_color;

    }

    for (auto &superColor : VQ.temp_superColors) {
        int Q = superColor.first;
        superColor.second.erase(curr_kmer_color);
        if (superColor.second.empty()) {
            continue;
        }

        uint64_t super_color_id = VQ.create_super_color(superColor.second);
        bool super_color_exist = (VQ.superColors[Q].find(super_color_id) != VQ.superColors[Q].end());

        if (super_color_exist) {
            VQ.superColorsCount[Q][super_color_id]++;
        } else {
            VQ.superColors[Q][super_color_id] = VQ.temp_superColors[Q];
            VQ.superColorsCount[Q][super_color_id] = 1;
        }

    }

    //     // Exporting superColors

    // for (auto &superColor : VQ.superColors) {
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

    // for (auto &superColor : VQ.superColorsCount) {
    //     int Q = superColor.first;
    //     cout << "Q" << Q << endl;
    //     for (auto const &color : superColor.second) {
    //         uint64_t supercolor_id = color.first;
    //         cout << supercolor_id << ": " << color.second;
    //         cout << endl;
    //     }
    //     cout << endl;
    // }


    return EXIT_SUCCESS;
}
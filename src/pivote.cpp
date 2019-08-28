#include "pivote.hpp"


void pivote(string index_prefex, set<int> allQs) {

    MAP2 edges;
    uint32_t seq1, seq2;
    uint16_t min_kmers, shared;
    uint16_t Q;
    string temp;

    string pairwise_file = index_prefex + "_kCluster.tsv";

    ifstream infile(pairwise_file);

    std::getline(infile, temp);

    while (infile >> seq1 >> seq2 >> min_kmers >> Q >> shared) {
        edges[{{seq1, seq2}, min_kmers}][Q] = shared;
        // edges[{{seq1,seq2},min_kmers}].emplace_back(Q, shared);
        // if(shared > max_shared) max_shared = shared;
        // cout << seq1 << "," << seq2 << "," << min_kmers << "," << Q << "," << shared << endl;
    }

//    cerr << "Size: " << edges.size() << endl;

    // EXPORTING
    string output_file = index_prefex + "_pivoted.tsv";
    ofstream pivoted;
    pivoted.open(output_file);

    string delimiter = "";

    for (auto const &edge : edges) {
        pivoted << "seq1\tseq2\tmin_kmers\t";

        for (auto const &Q : allQs) {
            pivoted << delimiter << "Q" << (int) Q;
            delimiter = "\t";
        }
        pivoted << "\n";

        break;
    }

    for (auto const &edge : edges) {

        pivoted << edge.first.first.first << '\t' << edge.first.first.second << '\t' << edge.first.second << '\t';
        delimiter = "";
        int i = 0;


        for (auto const &Q : allQs) {
            pivoted << delimiter;
            auto it = edge.second.find(Q);
            if (it != edge.second.end()) pivoted << it->second;
            else pivoted << 0;
            delimiter = "\t";
        }

        pivoted << endl;

    }

    pivoted.close();

}

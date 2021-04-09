#include "kSpider.hpp"

int main(int argc, char **argv) {

    if (argc != 2) {
        cerr << "must provide an index prefix\nrun: ./kSpider_pairwise <idx_prefix>" << endl;
        exit(1);
    }

    string index_prefix = argv[1];

    kSpider::pairwise(index_prefix);

}
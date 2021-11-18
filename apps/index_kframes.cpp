#include "kSpider.hpp"

int main(int argc, char **argv) {

    if (argc != 2) {
        cerr << "must provide an kframes directory\nrun: ./index_kframes <kfs_dir>" << endl;
        exit(1);
    }

    string kfs_dir = argv[1];

    kSpider::index_datasets(kfs_dir);

}
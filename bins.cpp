#include "kSpider.hpp"

inline uint64_t to_uint64_t(std::string const& value) {
    uint64_t result = 0;
    char const* p = value.c_str();
    char const* q = p + value.size();
    while (p < q) {
        result *= 10;
        result += *(p++) - '0';
    }
    return result;
}

int main(int argc, char** argv) {
    if(argc < 6){
        cout << "args: <bins_dir> <kSize> <output_prefix> <initial_reserve_size> <legend_reserve>\n";
        exit(1);
    }
    string bins_dir = argv[1];
    int kSize = stoi(argv[2]);
    string output_prefix = argv[3];
    uint64_t reserve_size = to_uint64_t(argv[4]);
    uint64_t legend_reserve = to_uint64_t(argv[5]);

    kSpider::bins_indexing(bins_dir, kSize, output_prefix, reserve_size, legend_reserve);
}
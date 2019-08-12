#include <iostream>
#include <vector>
#include <utility>
#include <stdint.h>

class Combo{

public:
    Combo(){}
    std::vector<std::pair<uint32_t , uint32_t>> combs;
    void combinations(int n);

private:
    int *arr = new int[2];
    int r = 2;
    void comb(int n, int r, int *arr);

};

#include "combinations.hpp"

void Combo::combinations(int n){
    this->combs.clear();
    this->comb(n, this->r, this->arr);
}

void Combo::comb(int n, int r, int *arr) {

    for (int i = n; i >= r; i--) {
        // choose the first element
        arr[r - 1] = i;
        if (r > 1) { // if still needs to choose
            // recursive into smaller problem
            comb(i - 1, r - 1, arr);

        } else {
            this->combs.emplace_back(std::make_pair(arr[0] - 1, arr[1] - 1));
        }
    }
}
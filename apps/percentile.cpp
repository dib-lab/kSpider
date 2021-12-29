#include <iostream>
#include <vector>
#include <cmath>
#include <bits/stdc++.h>

using namespace std;

int percentile(vector<uint64_t> & vec, uint64_t perc){
        sort(vec.begin(), vec.end());
        return vec[ceil((vec.size() * perc / 100))];
}

int main(int argc, char** argv) {
    vector<uint64_t> vec; 
    vector<uint64_t> a = { 1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,10,20,1000 };
    for(auto x: a) vec.push_back(x);
    cout << percentile(vec, stoi(argv[1])) << endl;
}
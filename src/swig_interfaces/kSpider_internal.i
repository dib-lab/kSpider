%module kSpider_internal

%{
#include "kSpider.hpp"
%}

using namespace std;
%include std_string.i

namespace kSpider{
    void pairwise(string index_prefix);
};
#pragma once
#include "types/integers.h"
#include "vector"
using namespace std;
template <integral T>
vector<T> all_divisor_naive(T a){
    vector<T> s;
    T i;
    for(i = 1; i * i < a; ++i) if(a % i == 0) s.push_back(i);
    int n = s.size();
    if(i * i == a) s.push_back(i);
    for(int j = n-1; j >= 0; --j) s.push_back(a/s[j]);
    return s;
}
#pragma once
#include <vector>
#include <iostream>
#include <cassert>
using namespace std;
template <typename T>
ostream& operator<<(ostream& f, const vector<T>& v){
    for(auto &i: v) f << i << ' ';
    return f;
}

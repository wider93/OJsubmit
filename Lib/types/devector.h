#pragma once
#include <vector>

using namespace std;
template <typename T>
struct DEdeVector {
    vector<T> pos, neg;
    explicit DEdeVector(int n, int m=0, const T& v = 0): pos(n, v), neg(m, v){}
    int sz() const { return psz()+nsz(); }
    int psz() const { return pos.size(); }
    int nsz() const { return neg.size(); }
    T& operator[](int n) {
        if (n >= 0) return pos[n];
        return neg[~n];
    }
    const T& operator[](int n) const {
        if (n >= 0) return pos[n];
        return neg[~n];
    }
    T get(int n) const {
        if(-nsz() <= n && n < psz()) return (*this)[n];
        return T(0);
    }
    void resize(int n, const T& val = 0) {
        pos.resize(n, val); neg.resize(n, val);
    }
    void reduce(){
        while(!pos.empty() && !pos.back()) pos.pop_back();
        while(!neg.empty() && !neg.back()) neg.pop_back();
    }
};

#pragma once
#include <iostream>
#include <array>
#include <algorithm>
#include <utility>
#include <initializer_list>
using namespace std;
template <typename T, int N>
struct avector{
    int sz;
    array<T, N> v;
    constexpr avector(int n=0, T val=T{}): sz(n), v{} {
        fill(v.begin(), v.begin() + sz, val);
    }
    constexpr avector(initializer_list<T> l):sz(l.size()), v{} {
        move(l.begin(), l.end(), v.begin());
    }
    using iterator = typename array<T, N>::iterator;
    constexpr unsigned long long size() const { return sz;}
    constexpr long long ssize() const { return sz;}
    constexpr void resize(int n) { sz = n;}
    constexpr auto begin() {return v.begin();}
    constexpr auto begin() const {return v.begin();}
    constexpr auto end() { return v.begin()+sz;}
    constexpr auto end() const { return v.begin()+sz;}
    constexpr auto cbegin() const { return v.cbegin();}
    constexpr auto cend() const { return v.cbegin()+sz;}
    constexpr auto& back() {return *(v.begin() + sz - 1);}
    constexpr const auto& back() const { return *(v.begin() + sz - 1);}
    constexpr bool empty() const { return sz == 0;}
    constexpr const T& operator[] (int i) const { return v[i];}
    constexpr T& operator[](int i) { return v[i];}
    constexpr void push_back(const T& a){ v[sz++] = a;}
    template<typename... Args>
    constexpr void emplace_back(Args&&... args ){
        v[sz++] = T(forward<Args>(args)...);
    }
    constexpr weak_ordering	operator <=> (const avector& o) const {
        for(int i = 0, msz = min(sz, o.sz); i < msz; ++i) if(v[i] != o.v[i]) return v[i] <=> o.v[i];
        return sz <=> o.sz;
    }
    constexpr bool operator!= (const avector& o) const {
        if (sz != o.sz) return true;
        for(int i = 0; i < sz; ++i) if(v[i] != o.v[i]) return true;
        return false;
    }
    constexpr bool operator== (const avector& o) const { return !((*this) != o);}
    constexpr void pop_back(){--sz;}
};

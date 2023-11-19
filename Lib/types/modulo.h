#pragma once
#include "types/integers.h"
#include <iostream>
#include <vector>
#include <traits.h>

using namespace std;
template <i64 pv>
struct Z_n{
    using T = conditional_t<(pv > 2147483647), i64, i32>;
    using Tw = conditional_t<(pv > 46340), conditional_t<(pv > 2147483647), i128, i64>, i32>;
    constexpr static T mod = pv;
public:
    T v;
    constexpr static i64 getmod() {return mod;}
    constexpr Z_n(const T& a = 0) noexcept: v(a){}
    explicit constexpr operator bool() const{ return v != 0; }
    template<integral I>explicit constexpr operator I() const{ return v; }
    friend istream& operator >> (istream &in, Z_n<mod> &a) { in >> a.v; a.v %= mod; if(a.v < 0) a.v += mod; return in;}
    friend ostream& operator << (ostream &out, const Z_n<mod> &a) { return out << a.v; }
    constexpr Z_n  operator+ (const Z_n & o) const{ T tmp = v - mod + o.v; if(tmp < 0) tmp += mod; return tmp; }
    constexpr Z_n& operator+=(const Z_n & o){ v += o.v - mod; if(v < 0) v += mod; return *this; }
    constexpr Z_n& operator++(){return *this += 1;}
    constexpr Z_n  operator++(int) {Z_n ans = *this; *this += 1; return ans;}
    constexpr Z_n  operator- (const Z_n & o) const{T tmp = v - o.v; if(tmp < 0) tmp += mod; return tmp; }
    constexpr Z_n& operator-=(const Z_n & o){ v -= o.v; if(v < 0) v += mod; return *this; }
    constexpr Z_n& operator--(){return *this -= 1;}
    constexpr Z_n  operator--(int) { Z_n ans = *this; *this -= 1; return ans; }
    constexpr Z_n  operator* (const Z_n & o) const{ return (Tw)v * o.v % mod; }
    constexpr Z_n& operator*=(const Z_n & o){ v = (Tw)v * o.v % mod; return *this; }
    constexpr Z_n  operator- () const{ return v == 0 ? 0 : mod - v; }
    constexpr Z_n  operator/ (const Z_n & o) const{ return *this * o.inv(); }
    constexpr Z_n& operator/=(const Z_n & o){ return *this = *this * o.inv(); }
    constexpr bool operator==(const Z_n & o) const{ return v == o.v; }
    constexpr bool operator!=(const Z_n & o) const{ return v != o.v; }
    constexpr Z_n& iaddmul(const Z_n& a, const Z_n& b) {v = (v + (long long)a.v * b.v) % mod; return *this;}
    constexpr Z_n addmul(const Z_n& a, const Z_n& b) const {return (v + (long long)a.v * b.v) % mod;}
    constexpr Z_n pow(u64 exp) const{ return pow(*this, exp); }
    constexpr static Z_n pow(Z_n base, u64 exp){
        Z_n ans = 1;
        while(exp ){
            if(exp & 1) ans *= base;
            base *= base;
            exp >>= 1;
        } return ans;
    }
    constexpr Z_n inv() const{ return inv(v); }
    constexpr static T inv(T x) {
        T a = egcd(x, mod);
        return (a < 0) ? a + mod : a;
    }
};
namespace std{template<i64 p>struct hash<Z_n<p>>{u64 operator()(const Z_n<p> &a)const {return a.v;}};}
template <typename T>
concept is_modular = Number<T> && requires(T a) {
    {Z_n<T::mod>(a)} -> same_as<T>;
};

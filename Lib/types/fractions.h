#pragma once
#include "traits.h"
#include <numeric>
using namespace std;
template <integral T>
struct Fraction{
    T m, n;
    constexpr Fraction(T m = 0, T n = 1):m(m), n(n){}
    constexpr Fraction& reduce(bool sign= false){
        if(sign && n < 0){m = -m, n = -n;}
        auto g = gcd(m, n); // g = gcd(abs(m), abs(n))
        m /= g; n /= g;
        return *this;
    }
    constexpr Fraction operator+(const Fraction& o) const{ return {m * o.n + o.m * n, n * o.n}; }
    constexpr Fraction operator-(const Fraction& o) const{ return {m * o.n - o.m * n, n * o.n}; }
    constexpr Fraction operator*(const Fraction& o) const{ return {m * o.m, n * o.n}; }
    constexpr Fraction operator/(const Fraction& o) const{ return {m * o.n, n * o.m}; }
    constexpr Fraction& operator+=(const Fraction& o){ return *this = *this + o; }
    constexpr Fraction& operator-=(const Fraction& o){ return *this = *this - o; }
    constexpr Fraction& operator*=(const Fraction& o){ return *this = *this * o; }
    constexpr Fraction& operator/=(const Fraction& o){ return *this = *this / o; }
    constexpr bool operator<(const Fraction &o) const {return m * o.n < o.m * n;}
};

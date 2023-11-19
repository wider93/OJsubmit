#pragma once
#include <iostream>
#include <numeric>
using namespace std;
template <integral T, bool autoReduction = false >
struct Fraction{
    T m, n;
    constexpr Fraction(T m = 0, T n = 1):m(m), n(n){if constexpr (autoReduction) reduce();}
    constexpr Fraction& reduce(bool sign= false){
        if(sign && n < 0){m = -m, n = -n;}
        auto g = gcd(m, n); // g = gcd(abs(m), abs(n))
        m /= g; n /= g;
        return *this;
    }
    constexpr Fraction abs() const{return {std::abs(m), n};}
    constexpr Fraction operator+(const Fraction& o) const{ return {m * o.n + o.m * n, n * o.n}; }
    constexpr Fraction operator-(const Fraction& o) const{ return {m * o.n - o.m * n, n * o.n}; }
    constexpr Fraction operator*(const Fraction& o) const{ return {m * o.m, n * o.n}; }
    constexpr Fraction operator/(const Fraction& o) const{ return {m * o.n, n * o.m}; }
    constexpr Fraction& operator+=(const Fraction& o){ return *this = *this + o; }
    constexpr Fraction& operator-=(const Fraction& o){ return *this = *this - o; }
    constexpr Fraction& operator*=(const Fraction& o){ return *this = *this * o; }
    constexpr Fraction& operator/=(const Fraction& o){ return *this = *this / o; }
    constexpr Fraction operator-()const { return {-m, n}; }
    constexpr auto operator<=>(const Fraction &o) const {return m * o.n <=> o.m * n;}
    constexpr auto operator==(const Fraction &o) const {return m * o.n == o.m * n;}
    friend ostream& operator<<(ostream& f, const Fraction& o){
        return f << o.m << '/' << o.n;
    }
    template<floating_point S>
    explicit operator S() const {return S(m)/n;}
};

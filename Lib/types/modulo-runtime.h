#pragma once
#include <types/integers.h>
#include <basic/utility.h>
#include <bit>
using namespace std;

template <int inputLength, int modLength = inputLength>
requires requires{ inputLength >= modLength; inputLength <= 64;}
struct Barret{
    using T = conditional_t<inputLength < 31, i32, i64>;
    using Tw = conditional_t<modLength + inputLength <= 63, i64, i128>;
    T mod, m;
    constexpr Barret(T n=2): mod(n), m(((Tw)1 << inputLength) / n){
        assert((mod >> m) == 0 && mod > 0);
    }
    friend istream& operator>>(istream& f, Barret& o){ T a; f >> a; o = Barret(a); return f;}
    friend T operator/(Tw a, const Barret& d) {
        Tw q = a * d.m >> inputLength;
        a -= q * d.mod;
        if(a >= d.mod) ++q;
        return q;
    }
    friend T operator%(Tw a, const Barret& d) {
        Tw q = a * d.m >> inputLength;
        a -= q * d.mod;
        if(a >= d.mod) a -= d.mod;
        return a;
    }
    friend T& operator/=(T& a, const Barret& d) { return a = a / d; }
    friend T& operator%=(T& a, const Barret& d) { return a = a % d; }
};
template <int k, int l=2*k, int index=0> // "index" exists solely to fulfill the requirement of keeping multiple modulos.
struct RZ_n{
    using B = Barret<l, k>;
    using T = typename B::T;
    using Tw = wider_t<T>;
    inline static B mod{3};
    T v;
    constexpr RZ_n(const T& a = 0) noexcept: v(a){}
    static i64 getmod() {return mod.mod;}
    static void setmod(typename B::T a) {mod = B(a);}
    explicit constexpr operator bool() const{ return v != 0; }
    template<integral I>explicit constexpr operator I() const{ return v; }
    friend istream& operator>> (istream &in, RZ_n &a) { in >> a.v; a.v %= mod.mod; if(a.v < 0) a.v += mod.mod; return in;}
    friend ostream& operator<< (ostream &out, const RZ_n &a) { return out << a.v; }
    constexpr RZ_n  operator+ (const RZ_n & o) const{ T tmp = v - mod.mod + o.v; if(tmp < 0) tmp += mod.mod; return tmp; }
    constexpr RZ_n& operator+=(const RZ_n & o){ v += o.v - mod.mod; if(v < 0) v += mod.mod; return *this; }
    constexpr RZ_n& operator++(){return *this += 1;}
    constexpr RZ_n  operator++(int) {RZ_n ans = *this; *this += 1; return ans;}
    constexpr RZ_n  operator- (const RZ_n & o) const{T tmp = v - o.v; if(tmp < 0) tmp += mod.mod; return tmp; }
    constexpr RZ_n& operator-=(const RZ_n & o){ v -= o.v; if(v < 0) v += mod.mod; return *this; }
    constexpr RZ_n& operator--(){return *this -= 1;}
    constexpr RZ_n  operator--(int) { RZ_n ans = *this; *this -= 1; return ans; }
    constexpr RZ_n  operator* (const RZ_n & o) const{ return (Tw)v * o.v % mod; }
    constexpr RZ_n& operator*=(const RZ_n & o){ v = (Tw)v * o.v % mod; return *this; }
    constexpr RZ_n  operator- () const{ return v == 0 ? 0 : mod.mod - v; }
    constexpr RZ_n  operator/ (const RZ_n & o) const{ return *this * o.inv(); }
    constexpr RZ_n& operator/=(const RZ_n & o){ return *this = *this * o.inv(); }
    constexpr bool operator==(const RZ_n & o) const{ return v == o.v; }
    constexpr bool operator!=(const RZ_n & o) const{ return v != o.v; }
    constexpr RZ_n pow(unsigned long long exp) const{ return pow(*this, exp); }
    constexpr static RZ_n pow(RZ_n base, unsigned long long exp){
        RZ_n ans = 1;
        while(exp){
            if(exp & 1) ans *= base;
            base *= base;
            exp >>= 1;
        } return ans;
    }
    constexpr RZ_n inv() const{ return inv(v); }
    constexpr static T inv(T x) {
        T a = egcd(x, mod.mod);
        return (a < 0) ? a + mod.mod : a;
    }
};
template <int maxmod>
using RZ_from_mod = RZ_n<bit_width((u64)maxmod)>;
template <typename T>
concept is_runtime_modular = Number<T> && requires(T a) {
    typename T::B;
//    {RZ_n<T::B::k, T::index>(a)} -> same_as<T>;
};
#pragma once
#include <iostream>
#include <utility>

using namespace std;

using i32 = int;
using i64 = long long;
using i128 = __int128;
using u32 = unsigned int;
using u64 = unsigned long long;
using u128 = unsigned __int128;

template<integral T> struct wider { using type = T; };
template<> struct wider<i32> { using type = i64; static constexpr i32 width = 32;};
template<> struct wider<i64> { using type = i128;static constexpr i32 width = 64;};
template<> struct wider<u32> { using type = u64; static constexpr i32 width = 32;};
template<> struct wider<u64> { using type = u128;static constexpr i32 width = 64;};

template <integral T> using wider_t = typename wider<T>::type;
template <integral T>
inline constexpr T mulmod(T a, T b, T c){return wider_t<T>(a) * b % c;}
template <integral T>
inline constexpr T pow(T base, long long exp, T p){
    T ans = 1;
    while(exp){
        if(exp & 1) ans = mulmod(ans, base, p);
        base = mulmod(base, base, p);
        exp >>= 1;
    } return ans;
}
template <integral T, integral S>
inline constexpr auto pow(T base, long long exp, S p){
    return pow<common_type_t<T,S>>(base, exp, p);
}
template <signed_integral T>
constexpr T divint(T a, T b){ // b > 0
    return (a >= T(0)) ? a/b : ~((~a)/b);
}
template <signed_integral T>
constexpr T egcd(T x, T y) {
    T a = 1, b = 0;
    if(x < 0) x = -x, a = -a;
    while(x){
        T d = y / x;
        y = exchange(x, y - d * x);
        b = exchange(a, b - d * a);
    }  return b;
}
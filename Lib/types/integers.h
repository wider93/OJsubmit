#pragma once
#include <iostream>
#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>
#include <cassert>
#include <bit>

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
constexpr T divint(T a, T b){
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
template <unsigned_integral T>
struct OddMont{
    using S = wider_t<T>;
    const T n, ni, r2, r1, r3;
    static constexpr T inv_word(T a){ // inverse of 1 << wordsize. a must be odd.
        T ans = a;
        for(int i = 0; i < 5; ++i) ans *= (T)2 - ans * a;
        return ans;
    }
    static constexpr T redc_given(S a, T n, T ni) {
        T m = (a & (T)-1)* ni;
        if(a >= S(-m)*n){
            return (a - (S(-m)) * n) >> wider<T>::width;
        } return (a + (S)m*n) >> wider<T>::width;
    }
    constexpr OddMont(T a):n(a), ni(-inv_word(a)), r2(-(S)a % a), r1(redc_given(r2, n, ni)), r3(redc_given(r2*r2, n, ni)) {}
    struct M{T v;}; // using seperate type to distinguish montgomery forms and normal integers
    constexpr M  redc(S a) const { return {redc_given(a, n, ni)};}
    constexpr M   add(M  a, M b) const{ if(a.v >= n - b.v) return {a.v+b.v-n}; return {a.v+b.v};}
    constexpr M& iadd(M &a, M b) const{ a.v += b.v-n; if(b.v > a.v) a.v += n; return a;}
    constexpr M   sub(M  a, M b) const{ if(a.v > b.v) return {n - b.v+a.v}; return {a.v - b.v};}
    constexpr M& isub(M &a, M b) const{ if(a.v > b.v) a.v=n - b.v+a.v; else a.v -= b.v; return a;}
    constexpr M   mul(M  a, M b) const{ return redc(S(a.v) * b.v);}
    constexpr M& imul(M &a, M b) const{ return a = redc(S(a.v) * b.v);}
    constexpr M toMont(T a) const {return redc(S(r2) * a);}
    constexpr M   add(M  a, T b) const{ return  add(a, toMont(b));}
    constexpr M& iadd(M &a, T b) const{ return iadd(a, toMont(b));}
    constexpr M   sub(M  a, T b) const{ return  sub(a, toMont(b));}
    constexpr M& isub(M &a, T b) const{ return isub(a, toMont(b));}
    constexpr M   mul(M  a, T b) const{ return  mul(a, toMont(b));}
    constexpr M& imul(M &a, T b) const{ return imul(a, toMont(b));}
    constexpr T powmul(M base, long long exp, T v) const{
        M ans = {v};
        while(exp){
            if(exp & 1) imul(ans, base);
            imul(base, base);
            exp >>= 1;
        } return ans.v;
    }
    constexpr T pow(M base, long long exp) const{
        return powmul(base, exp, 1);
    }
    constexpr T pow(T base, long long exp) const{
        return pow(toMont(base), exp);
    }
};

constexpr bool is_prime_naive(i64 k){
    if(k < 2) return false;
    if(k % 2 == 0) return k == 2;
    for(i64 q = 3; q*q <= k; q += 2) if(k%q==0) return false;
    return true;
}
namespace Details{
    constexpr int bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    constexpr pair<u64, i32> threshold_array[] = {{3825123056546413051ull, 9}, {341550071728321ull, 7}, {3474749660383ull, 6}, {2152302898747ull, 5}, {3215031751ull, 4}};
    constexpr i32 thres(u64 k){
        i32 t = 12;
        for(const auto &[i, j]: threshold_array){
            if(k < i) t = j;
            else break;
        }return t;
    }
    constexpr bool miller_rabin_base(const OddMont<u64>& k, i32 s, u64 odd, i32 base){
        //inputLength - 1 = odd << s
        typename OddMont<u64>::M x = {k.powmul(k.toMont(base), odd, k.r1)};
        for(int i = 0; i < s && x.v != k.r1; ++i){
            auto y = k.mul(x, x);
            if(y.v == k.r1 && x.v != k.n - k.r1) return false;
            x = y;
        } return x.v == k.r1;
    }
    constexpr bool miller_rabin(u64 k){
        if(k % 2 == 0) return k == 2;
        i32 s = 0; u64 odd = k-1;
        while((odd & 1) == 0) {odd >>= 1; ++s;}
        auto mod = OddMont<u64>(k);
        for(i32 i = 0; i < thres(k); ++i){
            if(bases[i] == k) return true;
            if(!miller_rabin_base(mod, s, odd, bases[i])) return false;
        } return true;
    }
}
constexpr bool is_prime(i64 k){
    if(k < 1000) return is_prime_naive(k);
    return Details::miller_rabin(k);
}

template <unsigned_integral T>
constexpr T pollard_rho_type(const T r, T init=1) {
    using S = wider_t<T>;
    const auto rm = OddMont<T>(r);
    using M = typename decltype(rm)::M;
    T x0 = 2, k = r, itr;
    const u32 step = 1u << 10;
    while (k == r) {
        M y{x0++}, x=y;
        for(k = 1, itr = step; k == 1 && itr < (1u << 21); itr <<= 1){
            M g{1};
            x = y;
            for(T s = 0; s < itr && k == 1; s += step){
                for(T i = 0; i < step; ++i){
                    y = rm.redc((S)(y.v)*y.v+init);
                    g = rm.imul(g, max(x.v, y.v) - min(x.v, y.v));
                } k = gcd(g.v, r);
                if (k == r){
                    k = 1; auto py = x;
                    for(T i = 0; i < step && k == 1; ++i){
                        py = rm.redc((S)(py.v)*py.v+init);
                        k = gcd(r, max(x.v, py.v) - min(x.v, py.v));
                    }if(k == 1) k = r;
                }
            }
        }
    } return k;
}
constexpr i64 pollard_rho(u64 r){
    if(r < (1ull << 32)) return pollard_rho_type<u32>(r);
    return pollard_rho_type<u64>(r);
}
template <typename T=u64, typename Container = vector<T>>
constexpr Container factorize(T n){
    Container v;
    for(const auto &p: {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}){
        while(n % p == 0){
            v.push_back(p);
            n /= p;
        }
    }if(n == 1) return v;
    vector<pair<T,i32>> q;
    q.emplace_back(n, 1);
    while(!q.empty()){
        auto [t, num] = q.back();
        q.pop_back();
        if(t < 1681 or Details::miller_rabin(t)) for(;num;--num)v.push_back(t);
        else{
            {   i64 s = sqrtl(t)+.5;
                if(s * s == t) {
                    q.emplace_back(s, num<<1); continue;
                }
            }
            auto a = pollard_rho(t);
            {q.emplace_back(a, num); q.emplace_back(t/a, num);}
        }
    }
    sort(begin(v), end(v));
    return v;
}
template <integral T>
inline constexpr T primitive_root(T prime){
    if(prime == 2) return 1;
    auto factor = factorize(prime-1);
    for(T a = 2; a < prime; ++a){
        bool f = true;
        for(auto &p: factor){
            if(pow(a, (prime - 1) / p, prime) == (T)1){
                f = false; break;
            }
        } if(f) return a;
    } assert(false);
}
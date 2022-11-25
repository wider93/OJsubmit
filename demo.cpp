#include <iostream>
#include <numeric>
#include <vector>
#include "sub.h"
#include <algorithm>
#include <cmath>

using namespace std;

using i32 = int;
using i64 = long long;
using i128 = __int128;
using u32 = unsigned int;
using u64 = unsigned long long;
using u128 = unsigned __int128;

template<integral T> struct wider { using type = T; };
template<> struct wider<i32> { using type = i64; };
template<> struct wider<i64> { using type = i128; };
template<> struct wider<u32> { using type = u64; };
template<> struct wider<u64> { using type = u128; };

template <typename T> using wider_t = typename wider<T>::type;
template <typename T>
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

constexpr bool is_prime_naive(i64 k){
    if(k < 2) return false;
    if(k % 2 == 0) return k == 2;
    for(i64 q = 3; q*q <= k; q += 2) if(k%q==0) return false;
    return true;
}
namespace Details{
    constexpr int bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    constexpr pair<i64, i32> threshold_array[] = {{3825123056546413051ll, 9}, {341550071728321ll, 7}, {3474749660383ll, 6}, {2152302898747ll, 5}, {3215031751ll, 4}};
    constexpr i32 thres(i64 k){
        i32 t = 12;
        for(const auto &[i, j]: threshold_array){
            if(k < i) t = j;
            else break;
        }return t;
    }
    constexpr bool miller_rabin_base(i64 k, i32 s, i64 odd, i64 base){
        //k - 1 = odd << s
        i64 x = pow(base, odd, k);
        for(int i = 0; i < s && x != 1; ++i){
            i64 y = mulmod(x, x, k);
            if(y == 1 && x != k - 1) return false;
            x = y;
        } return x == 1;
    }
    constexpr bool miller_rabin(i64 k){
        if(k % 2 == 0) return k == 2;
        i32 s = 0; i64 odd = k-1;
        while((odd & 1) == 0) {odd >>= 1; ++s;}
        for(int i = 0; i < thres(k); ++i){
            if(bases[i] == k) return true;
            if(!miller_rabin_base(k, s, odd, bases[i])) return false;
        } return true;
    }
}
constexpr bool is_prime(i64 k){
    if(k < 1000) return is_prime_naive(k);
    return Details::miller_rabin(k);
}

constexpr i64 primitive_root(i64 p){
    i64 p_div[20]{}, p_num = 0; // 2 * 3 * ... * 29 > 2 ** 31
    i64 k = p-1;
    for(i64 i = 2; i * i <= k; ++i)if(k % i == 0){
        p_div[p_num++] = i;
        do{ k /= i; } while(k % i == 0);
    } if (k > 1) p_div[p_num++] = k;
    for(i64 a = 2; a < p; ++a){
        bool f = true;
        for(int i = 0; i < p_num; ++i) {
            i64 q = p_div[i];
            if(pow(a, (p-1)/q, p) == 1){f = false; break;}
        } if(f) return a;
    }return 1;
}

template <integral T>
inline T pollard_rho_type(const T r, T init) {
    using S = wider_t<T>;
    while (1) {
        T _x = 2, _y = 2;
        T k = 1;
        T itr = 2;
        _x = ((S)_x * _x + init) % r;
        _y = ((S)_y * _y + init) % r;
        _y = ((S)_y * _y + init) % r;
        k = gcd(abs(_x - _y), r);
        if(k != 1) {
            if (k == r) ++init;
            else return k;
        }
        _x = ((S)_x * _x + init) % r;
        _y = ((S)_y * _y + init) % r;
        _y = ((S)_y * _y + init) % r;
        k = gcd(abs(_x - _y), r);
        if(k != 1) {
            if (k == r) ++init;
            else return k;
        }
        while (k == 1) {
            T g = 1;
            for(T i = 0; i < itr; ++i){
                _x = ((S)_x * _x + init) % r;
                _y = ((S)_y * _y + init) % r;
                _y = ((S)_y * _y + init) % r;
                g = (S)g * abs(_x - _y) % r;
            } k = gcd(g, r);
            itr <<= 1;
        }
        if (k == r) ++init;
        else return k;
    }
}
i64 pollard_rho(i64 r, i64 init){
    if(r < (1ll << 31) && init < (1ll << 31)) return pollard_rho_type<i32>(r, init);
    return pollard_rho_type(r, init);
}
template <typename T=i64, typename Container = vector<T>>
Container factorize(i64 n){
    Container v;
    for(const auto &p: {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}){
        while(n % p == 0){
            v.push_back(p);
            n /= p;
        }
    }if(n == 1) return v;
    vector<pair<T,i32>> q;
    q.emplace_back(n, 1);
    i32 init = 0;
    while(!q.empty()){
        auto [t, num] = q.back();
        q.pop_back();
        if(t < 1681 or Details::miller_rabin(t)) for(;num;--num)v.push_back(t);
        else{
            {
                i64 s = sqrt(t)+.5; 
                if(s * s == t) {
                    q.emplace_back(s, num<<1); continue;
                }
            }
            auto a = pollard_rho(t, ++init);
            for(;num;--num){q.emplace_back(a, 1); q.emplace_back(t/a, 1);}
        }
    }
    sort(begin(v), end(v));
    return v;
}



int main(){
    cin.tie(0) -> sync_with_stdio(0);
    i32 tc; cin >> tc;
    i32 ans = 0;
    for(i32 i = 0; i < tc; ++i){
        i64 f; cin >> f;
        auto t = factorize(f);
        cout << t.size() << ' ';
        for(const auto& i: t) cout << i << ' ';
        cout << '\n';
    }
}


#pragma once
#include "types/integers.h"
#include "cassert"
#include "vector"
#include "numeric"
#include "cmath"
#include "algorithm"
#include "bitset"
using namespace std;
template <unsigned_integral T>
struct OddMont{
    using S = wider_t<T>;
    const T n, ni, r2, r1;
    static constexpr T inv_word(T a){ // inverse of 1 << word_size. a must be odd.
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
    constexpr OddMont(T a):n(a), ni(-inv_word(a)), r2(-(S)a % a), r1(redc_given(r2, n, ni)) {}
    struct M{T v;}; // using separate type to distinguish montgomery forms and normal integers
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
    constexpr i32 bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    constexpr pair<u64, i32> threshold_array[] = {{3825123056546413051ull, 9}, {341550071728321ull, 7}, {3474749660383ull, 6}, {2152302898747ull, 5}, {3215031751ull, 4}};
    constexpr i32 threshold(u64 k){
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
        for(i32 i = 0; i < threshold(k); ++i){
            if((u32)bases[i] == k) return true;
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
template <typename T=u64, typename Container = vector<pair<T, i32>>>
constexpr Container factorize_concat(T n){
    Container ans;
    auto v = factorize(n); auto m = ssize(v);
    for(int i = 0; i < m;){
        int j = i+1;
        for(;j < m && v[j] == v[j-1]; ++j);
        ans.emplace_back(v[i], j-i);
        i = j;
    } return ans;
};
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
#define eratosthenes_macro(n) sieve[0] = sieve[1] = false; \
    for(int i = 2; i * i <= n; ++i) if(sieve[i]) for(int j = i * i; j <= n; j += i) sieve[j] = false; \
    return sieve;

vector<char> sieve_eratosthenes(int n){
    assert(n > 0);
    vector<char> sieve(n+1, true);
    eratosthenes_macro(n)
}
template <int N>
bitset<N> sieve_eratosthenes_bitset(){
    bitset<N> sieve;
    sieve.flip();
    eratosthenes_macro(N);
}
auto sieve_eratosthenes_factor(int n){
    assert(n > 0);
    vector<int> sieve(n+1);
    for(int i = 0; i <= n; ++i) sieve[i] = i;
    for(int i = 2; i * i <= n; ++i) if(sieve[i] == i) for(int j = i * i; j <= n; j += i) sieve[j] = i;
    return sieve;
}
constexpr int fs = 64 * 3 * 5 * 7 * 11;
void apply_on_odd_primes(int n, auto& init, auto f) {
	if(n <= fs){
		auto sieve = sieve_eratosthenes_bitset<fs>();
		for(int i = 3; i <= n; i += 2) if(sieve[i]) { f(init, i); }
		return;
	}
	int cut = 1+sqrt(n-1);
	int sieve_cut = (cut + fs - 1) / fs * fs;
	vector<int> odd_primes;{
		auto sieve = sieve_eratosthenes(sieve_cut);
		for(int i = 3; i <= sieve_cut; i += 2) if(sieve[i]) {
			f(init, i);
			odd_primes.push_back(i);
		}
	}
	int s = odd_primes.size();
	while ((long long)odd_primes[s] * odd_primes[s] > n) --s;
	vector<int> current_rem(s);
	for(int i = 4; i < s; ++i) {
		current_rem[i] = odd_primes[i] - 1 - sieve_cut % odd_primes[i];
		if(current_rem[i] % 2) current_rem[i] += odd_primes[i];
		current_rem[i] /= 2;
	}
	bitset<fs/2> base; {
		for (int i: {3, 5, 7, 11})
			for (int j = i>>1; j < fs/2; j += i)
				base[j] = true;
	}
	int now = sieve_cut;
	for( ;now + fs < n; now += fs){
		bitset<fs/2> small = base;
		for(int i = 4; i < s; ++i){
			int j = current_rem[i];
			for(; j < fs/2; j += odd_primes[i])
				small[j] = true;
			current_rem[i] = j - fs/2;
		}
		for(int i = 0; i < fs/2; i++)
			if(!small[i])
				f(init, now + 2*i+1);
	}{
		for(int i = 4; i < s; ++i){
			for(int j = current_rem[i]; j < fs/2; j += odd_primes[i])
				base[j] = true;
		}for(int i = 1; i <= n-now; i += 2)
			if(!base[i>>1])
				f(init, now + i);
	}
}
#pragma once
#include "poly/fft.h"
using namespace std;
template<Number T, MulType use_FFT = MulType::naive>
struct Poly{
    vector<T> a;
    Poly() : a() {}
    Poly(T a0) : a(1, a0) { normalize(); }
    Poly(vector<T> &v) : a(v) { normalize(); }
    Poly(vector<T> &&v) : a(move(v)) { normalize(); }
    template<bidirectional_iterator A>
    Poly(A b, A e) : a(b, e){}
    Poly(initializer_list<T> x):a(begin(x), end(x)){}
    constexpr static Poly monomial(int n, T a = 1){
        vector<T> v(n+1); v[n] = a; return v;
    }
    constexpr void normalize(){ while(!a.empty() && a.back() == T(0)) a.pop_back(); }
    int size() const { return a.size(); }
    int deg() const { return size() - 1; }
    void push_back(const T &x){ a.push_back(x); }
    T get(int idx) const { return 0 <= idx && idx < size() ? a[idx] : T(0); }
    const T& operator[](int idx) const { return a[idx]; }
    T& operator[](int idx){ return a[idx]; }

    constexpr T leading() const { return size() ? a.back() : T(0); }
    constexpr Poly reversed() const { return {a.rbegin(), a.rend()}; }
    constexpr Poly trim(int sz) const { return {a.begin(), a.begin() + min(sz, size())}; }
    constexpr Poly inv_naive(int n) const {
        assert(get(0) != 0);
        T inv(T(1) / get(0));
        vector<T> b(n); b[0] = inv;
        for(int i = 1; i < n; ++i){
            T ans=0;
            for(int j = 0; j < i; ++j) ans -= get(i-j) * b[j];
            b[i] = ans * inv;
        } return b;
    }
    constexpr Poly inv(int n) const {
        assert(a[0] != 0);
        Poly q = inv_naive(min(32, n));
        for(int i = 32; i < n; i *= 2){
            Poly p = (q * trim(i * 2)).trim(i * 2);
            p[0] -= T(2);
            q = (p.neg() * q).trim(i * 2);
        } return q.trim(n);
    }

    constexpr Poly& operator*=(const T &x){ for(auto &i: a) i *= x; normalize(); return *this; }
    constexpr Poly  operator* (const T& x) const { return Poly(*this) *= x; }
    constexpr Poly& operator/=(const T &x){ T y = T(1)/x; for(auto &i: a) i *= y; return *this; }
    constexpr Poly  operator/ (const T& x) const { return Poly(*this) /= x; }

    constexpr Poly  operator+ (const Poly &o) const { return Poly(*this) += o; }
    constexpr Poly& operator+=(const Poly &o){
        a.resize(max(size(), o.size()));
        for(int i = 0; i < o.size(); i++) a[i] += o.a[i];
        normalize();
        return *this;
    }
    constexpr Poly& neg() {for(auto &i: a) i = -i; return *this;}
    constexpr Poly  operator- () const { return Poly(*this).neg();}
    constexpr Poly  operator- (const Poly &o) const { return Poly(*this) -= o; }
    constexpr Poly& operator-=(const Poly &o){
        if(o.size() > size()) a.resize(o.size());
        for(int i = 0; i < o.size(); i++) a[i] -= o.a[i];
        normalize();
        return *this;
    }
    constexpr Poly operator* (const Poly &o) const{
        return multiply<T, use_FFT>(a, o.a);
    }
    constexpr Poly operator*= (const Poly &o){ return *this = *this * o;}

    constexpr Poly operator/ (const Poly &o) const{
        if constexpr(use_FFT != MulType::naive){
            if(size() >= 128 && o.deg() > 1) return div_dnc(o);
        } return div_naive(o);
    }
    constexpr Poly& operator/= (const Poly &o){ return *this = *this / o; }

    constexpr Poly operator% (const Poly &o) const {
        if(deg() < o.deg()) return *this;
        if constexpr(use_FFT != MulType::naive){
            if(size() >= 128 && o.deg() > 1){
                return *this - (*this / o) * o;
            }
        }
        return mod_naive(o);
    }
    constexpr Poly& operator%= (const Poly &o){ return *this = *this % o; }
    constexpr Poly  operator<< (int n) const {vector<T> v(n+size()); copy(begin(a), end(a), begin(v)+n); return v;}
    constexpr Poly  operator>> (int n) const {return {begin(a)+min(size(), n), end(a)};}
    constexpr Poly& operator<<=(int n) {return *this = *this << n;}
    constexpr Poly& operator>>=(int n) {return *this = *this >> n;}
    constexpr pair<Poly, Poly> divmod_naive (const Poly &o){
        if(size() < o.size()) return {Poly(), *this};
        vector<T> res(size() - o.size() + 1);
        while(size() >= o.size()){
            int now = size() - o.size();
            res[now] = a.back() / o.a.back();
            for(int i = 0; i + 1 < o.size(); ++i)
                a[i+now] -= res[now] * o.a[i];
            a.pop_back();
        } return {res, *this};
    }
    constexpr Poly div_naive (const Poly &o) const { auto me = *this; return me.divmod_naive(o).first; }
    constexpr Poly mod_naive (const Poly &o) const { auto me = *this; return me.divmod_naive(o).second; }
    constexpr Poly div_dnc(const Poly &o) const {
        if(size() < o.size()) return Poly();
        int sz = size() - o.size() + 1;
        Poly ra = reversed().trim(sz), rb = o.reversed().inv(sz);
        ra = (ra * rb).trim(sz);
        ra.a.resize(sz);
        return ra.reversed();
    }
    T eval(const T& x) const {
        T res{0};
        for(int i = deg(); i >= 0; i--) res = res * x + a[i];
        return res;
    }
    T composite(const Poly& o) const {
        Poly res;
        for(int i = deg(); i >= 0; i--) res = res * o + a[i];
        return res;
    }
    Poly derivative(int n) const {
        vector<T> res(a.begin()+1, min(a.begin() + 1 + n, a.end()));
        for(int i = 0; i < res.size(); ++i) res[i] *= T(i+1);
        return res;
    }
    Poly integral (int n) const{
        int m = min(n-1, size());
        if constexpr(is_modular<T> || is_runtime_modular<T>){
            vector<T> res(m + 1);
            if(m) res[1] = 1;
            for(int i = 2; i <= m; ++i) res[i] = -res[T::getmod() % i] * (T::getmod() / i);
            for(int i = 0; i < m; i++) res[i+1] *= a[i];
            return res;
        }
        vector<T> res(m + 1);
        for(int i = 0; i < m; i++) res[i+1] = a[i] / T(i+1);
        return res;
    }
    Poly cyclic_convolution(const Poly& o, int n) const {
        Poly ans = *this * o;
        for(int i = ans.size() - 1 ; i - n >= 0; --i){
            ans[i-n] += ans[i];
            ans[i] = 0;
        } return ans.trim(n);
    }
    Poly cyclic_convolution(const Poly& o) const {
        assert(size() == o.size());
        return cyclic_convolution(o, size());
    }
    Poly log(int n) const {  // assert(n > 0);
        return (derivative(n-1) * inv(n)).trim(n-1).integral(n);
    }
    Poly exp(int n) const {  //assert(n > 0 && a[0] == T(0));
        if(size() == 0) return T(1);
        Poly res(T(1));
        for(int i = 1; i < n; i <<= 1){
            auto t = trim(i << 1) - res.log(i << 1) + T(1);
            res = (res * t).trim(min(n, i << 1));
        } return res;
    }
    Poly pow(i64 n, i32 bound) const {
        if(n == 0) return {1};
        if(bound <= 0 || size() == 0) return {};
        int z = 0;
        while(z < size() && a[z] == T(0)) ++z;
        if constexpr (is_modular<T> || is_runtime_modular<T>) {
            assert(bound < T::getmod());
            if(bound <= (i128)n*z) return {};
            int k = bound - n*z, ns = n % T::getmod();
            Poly tmp = (*this) >> z;
            T o = tmp[0]; tmp /= o;
            return ((tmp.log(k) * T(ns)).exp(k) << (n * z)) * o.power(n);
        } assert(0);
        return {};
    }
    Poly taylor_shift(T c) const {
        if(!size()) return {};
        if(c == T(0)) return *this;
        Poly A(size(), 1), C(size(), 1);
        vector<T> fac(size()+1, 1);
        for(int i=0; i<size(); ++i){
            A[i] = a[i] * fac[i];
            fac[i+1] = fac[i] * (i+1);
        } T invfac = fac[size()].inv(), cp = c.power(size());
        for(int i=0; i<size(); ++i) {
            invfac *= size() - i;
            cp /= c;
            C[i] = cp * invfac;
        }
        Poly B; {
            Poly B2 = A * C;
            if(B2.size() >= size())
                B = {B2.a.begin()+deg(), B2.a.end()};
        }
        invfac = fac[B.size()].inv();
        for(int i = B.size(); i > 0; --i) {
            invfac *= i;
            B[i-1] *= invfac;
        }
        return B;
    }
    friend ostream& operator<<(ostream& f, const Poly &v){return f << v.a;}
};

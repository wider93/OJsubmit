#pragma once
#include <types/modulo.h>
#include <types/modulo-runtime.h>
#include <complex>

#define cpx complex<double>
inline namespace FFT{

template<is_modular F> inline constexpr F prim = primitive_root(F::mod);
template<> inline constexpr Z_n<998244353> prim<Z_n<998244353>> = 3; // merely for BOJ GCC 11 compilation

template <Number F> constexpr F root_of_unity(int n, bool inv){
    assert(false); return F{};
};
template <is_modular F> constexpr F root_of_unity(int n, bool inv){
    assert(F :: mod % n == 1);
    return inv ? prim<F>.pow((F::mod - 1) / n).inv() : prim<F>.pow((F::mod - 1) / n);
}
template <> cpx root_of_unity<cpx>(int n, bool inv){
    return polar<double>(1, (inv?-2:2) * numbers::pi / n);
}
template <Number F> vector<F> roots_of_unity(int n, int m, bool inv){
    assert(false); return {};
}
template <is_modular F> vector<F> roots_of_unity(int n, int m, bool inv){
    auto r = root_of_unity<F>(n, inv);
    vector<F> ans(m);
    ans[0] = 1; for(int i = 1; i < m; ++i) ans[i] = ans[i-1] * r;
    return ans;
}
template <> vector<cpx> roots_of_unity<cpx>(int n, int m, bool inv){
    vector<cpx> ans(m);
    ans[0] = 1;
    for(int i = 0; (1 << i) < m; ++i){
        auto r = root_of_unity<cpx>(n >> i, inv);
        for(int j = (1 << i); j < min(m, 1 << (i + 1)); ++j)
            ans[j] = ans[j - (1 << i)] * r;
    }
    return ans;
}
template <Number F> vector<F> roots_of_unity(int n, bool inv){/// returns ROU of length n/2. n must be power of 2
    return roots_of_unity<F>(n, n/2, inv);
}

template<Number F>
void fft(vector<F> &f, bool inv = 0) {
    int n = f.size(), j = 0;
    for (int i = 1; i < n; i++) {
        int bit = (n >> 1);
        while (j >= bit) j -= bit, bit >>= 1;
        j += bit;
        if (i < j) swap(f[i], f[j]);
    }
    vector<F> root = roots_of_unity<F>(n, inv);
    for (int i = 2; i <= n; i <<= 1) {
        int step = n / i;
        for (int j = 0; j < n; j += i) {
            for (int k = 0; k < i / 2; k++) {
                F v = f[j | k | i / 2] * root[step * k];
                f[j | k | i / 2] = f[j | k] - v; f[j | k] += v; 
            }
        }
    }
    if (inv) {
        for (int i = 0; i < n; i++) f[i] /= n;
    }
}
template<Number F>
void constexpr hadamard(vector<F> &f, bool inv = 0) {
    int n = f.size(), j = 0;
    for (int i = 2; i <= n; i <<= 1) {
        for (int j = 0; j < n; j += i) {
            for (int k = 0; k < i / 2; k++) {
                F v = f[j | k | i / 2];
                f[j | k | i / 2] = f[j | k] - v; f[j | k] += v; 
            }
        }
    }
    if (inv) for (int i = 0; i < n; i++) f[i] /= n;
}

// helper for fft related type conversion
enum class MulType{
    naive, ntt, floating
};
template <Number T>
constexpr vector<T> multiply_naive(const vector<T> &a, const vector<T>& b){
    vector<T> res(a.size() + b.size());
    for(int i = 0; i < a.size(); i++)
        for(int j = 0; j < b.size(); j++)
            res[i+j] += a[i] * b[j];
    return res;
}
template <is_modular T>
constexpr vector<T> multiply_ntt(const vector<T> &a, const vector<T>& b){
    int n = 1;
    while(n < a.size() + b.size()) n *= 2;
    vector<T> _a = a, _b = b;
    _a.resize(n); _b.resize(n);
    FFT::fft(_a); FFT::fft(_b);
    for (int i = 0; i < n; i++) _a[i] *= _b[i];
    FFT::fft(_a, 1);
    return _a;
}
template <Number T>
requires requires {requires is_modular<T> || is_runtime_modular<T>;}
constexpr vector<T> multiply_fft_cut(const vector<T> &a, const vector<T>& b){
    int n = 1;
    while(n < a.size() + b.size()) n *= 2;
    vector<cpx> _a(n), _b(n), r1(n), r2(n);
    for(int i=0; i<a.size(); i++){
        _a[i] = cpx((i32)a[i] >> 15, (i32)a[i] & 32767);
    }
    for(int i=0; i<b.size(); i++){
        _b[i] = cpx((i32)b[i] >> 15, (i32)b[i] & 32767);
    }
    FFT::fft(_a); FFT::fft(_b);
    for (int i = 0; i < n; i++) {
        int j = (i ? (n - i) : i);
        cpx ans1 = (_a[i] + conj(_a[j]));
        cpx ans2 = (_a[i] - conj(_a[j])) * cpx(0, -1);
        cpx ans3 = (_b[i] + conj(_b[j]));
        cpx ans4 = (_b[i] - conj(_b[j])) * cpx(0, -1);
        r1[i] = (ans3 + cpx(-ans4.imag(), ans4.real())) * ans1 / 4.;
        r2[i] = (ans3 + cpx(-ans4.imag(), ans4.real())) * ans2 / 4.;
    }
    FFT::fft(r1, true);
    FFT::fft(r2, true);
    vector<T> ret(n);
    for(int i=0; i<n; i++){
        i128 av = llround(r1[i].real());
        i128 bv = llround(r1[i].imag()) + llround(r2[i].real());
        i128 cv = llround(r2[i].imag());
        ret[i] = ((av << 30) + (bv << 15) + cv) % T::getmod();
    } return ret;
}
template <Number T, MulType S>
constexpr vector<T> multiply(const vector<T> &a, const vector<T>& o){
    if(a.size() + o.size() >= 256){
        if constexpr (S == MulType::ntt){
            return multiply_ntt(a, o);
        }else if constexpr (S == MulType::floating) {
            return multiply_fft_cut(a, o);
        }
    } return multiply_naive(a, o);
}

} //end fft

template<Number T, MulType use_FFT = MulType::naive>
struct Poly{
    vector<T> a;
    Poly() : a() {}
    Poly(T a0) : a(1, a0) { normalize(); }
    Poly(const vector<T> &v) : a(v) { normalize(); }

    template<typename A>
    Poly(A b, A e) : a(b, e){}
    Poly(vector<T> &&v) : a(move(v)) { normalize(); }
    constexpr static Poly monomial(int n, T a = 1){
        vector<T> v(n); v[n-1] = a; return v;
    }
    constexpr void normalize(){ while(!a.empty() && a.back() == T(0)) a.pop_back(); }

    int size() const { return a.size(); }
    int deg() const { return size() - 1; }
    void push_back(const T &x){ a.push_back(x); }
    T get(int idx) const { return 0 <= idx && idx < size() ? a[idx] : T(0); }
    const T& operator[](int idx) const { return a[idx]; }
    T& operator[](int idx){ return a[idx]; }

    constexpr Poly reversed() const { return vector<T>(a.rbegin(), a.rend()); }
    constexpr Poly trim(int sz) const { return vector<T>(a.begin(), a.begin() + min(sz, size())); }
    constexpr Poly inv(int n) const {
        assert(a[0] != 0);
        Poly q(T(1) / a[0]);
        for(int i = 1; i < n; i *= 2){
            Poly p = (Poly(2) - q * trim(i * 2)).trim(i * 2);
            q = (p * q).trim(i * 2);
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
    constexpr Poly  operator- () const { Poly res(*this); for(auto &i: res.a) i = -i; return res;}
    constexpr Poly  operator- (const Poly &o) const { return Poly(*this) -= o; }
    constexpr Poly& operator-=(const Poly &o){
        a.resize(max(size(), o.size()));
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
            if(size() >= 64 && o.deg() > 1) return div_dnc(o);
        } return div_naive(o);
    }
    constexpr Poly& operator/= (const Poly &o){ return *this = *this / o; }

    constexpr Poly operator% (const Poly &o) const {
        if(deg() < o.deg()) return *this;
        if constexpr(use_FFT != MulType::naive){
            if(size() >= 64 && o.deg() > 1){
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
    Poly derivative() const {
        vector<T> res(a.begin()+1, a.end());
        for(int i = 0; i + 1 < size(); ++i) res[i] *= T(i+1);
        return res;
    }
    Poly integral () const{
        if constexpr(is_modular<T> || is_runtime_modular<T>){
            vector<T> res(size() + 1);
            if(size()) res[1] = 1;
            for(int i = 2; i <= size(); ++i) res[i] = -res[T::mod % i] * (T::mod / i);
            for(int i = 0; i < size(); i++) res[i+1] *= a[i];
            return res;
        }
        vector<T> res(size() + 1);
        for(int i = 0; i < size(); i++) res[i+1] = a[i] / T(i+1);
        return res;
    }
    Poly cyclic_convolution(const Poly& o, int n) const {
        Poly ans = *this * o;
        for(int i = ans.size() - 1 ; i - n >= 0; --i){
            ans[i-n] += ans[i];
            ans[i] = 0;
        }
        return ans.trim(n);
    }
    Poly cyclic_convolution(const Poly& o) const {
        assert(size() == o.size());
        return cyclic_convolution(o, size());
    }
    Poly log(int n) const {  // assert(n > 0);
        return (derivative() * inv(n)).trim(n-1).integral();
    }
    Poly exp(int n) const {  //assert(n > 0 && a[0] == T(0));
        if(size() == 0) return T(1);
        Poly res(T(1));
        for(int i = 1; i < n; i <<= 1){
            auto t = trim(i << 1) + Poly(1) - res.log(i << 1);
            res = (res * t).trim(i << 1);
        }
        return res.trim(n);
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
            return ((tmp.log(k) * T(ns)).exp(k) << (n * z)) * o.pow(n);
        }
        assert(0);
        return {};
    }
    friend ostream& operator<<(ostream& f, const Poly &v){
        for(auto &i: v.a) f << i << ' ';
        return f;
    }
};

#undef cpx

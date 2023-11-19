#pragma once
#include <types/modulo.h>
#include <types/modulo-runtime.h>
#include <types/primality.h>
#include <complex>
using namespace std;

#define cpx complex<double>
inline namespace FFT{
template<is_modular F> inline constexpr F prim = primitive_root(F::mod);
template<> inline constexpr Z_n<998244353> prim<Z_n<998244353>> = 3; // merely for BOJ GCC 11 compilation

template <Number F> constexpr F root_of_unity(int n, bool inv){
    assert(false); return F{};
};
template <Number F> requires (is_modular<F> || is_runtime_modular<F>)
F root_of_unity(int n, bool inv){
    assert(F :: getmod() % n == 1);
    F p; if constexpr (is_modular<F>) p = prim<F>; else p = primitive_root(F::getmod());
    return inv ? p.pow((F::getmod() - 1) / n).inv() : p.pow((F::getmod() - 1) / n);
}
template <> cpx root_of_unity<cpx>(int n, bool inv){
    return polar<double>(1, (inv?-2:2) * numbers::pi / n);
}
template <Number F> vector<F> roots_of_unity(int n, int m, bool inv){
    assert(false); return {};
}
template <Number F> requires (is_modular<F> || is_runtime_modular<F>)
vector<F> roots_of_unity(int n, int m, bool inv){
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
    int n = ssize(f), j = 0;
    for (int i = 1; i < n; i++) {
        int bit = (n >> 1);
        while (j >= bit) j -= bit, bit >>= 1;
        j += bit;
        if (i < j) swap(f[i], f[j]);
    }
    const vector<F> rou = roots_of_unity<F>(f.size(), inv);
    for (int i = 2; i <= n; i <<= 1) {
        int step = n / i;
        for (int j = 0; j < n; j += i) {
            for (int k = 0; k < i / 2; k++) {
                F v = f[j | k | i / 2] * rou[step * k];
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
    int n = ssize(f);
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

/** helper for fft related type conversion */
enum class MulType{
    naive, ntt, floating
};
template <Number T>
constexpr vector<T> multiply_naive(const vector<T> &a, const vector<T>& b){
    vector<T> res(max<int>(0, ssize(a) + ssize(b) - 1));
    for(int i = 0; i < ssize(a); i++)
        for(int j = 0; j < ssize(b); j++)
            res[i+j] += a[i] * b[j];
    return res;
}
template <Number T> requires requires {requires is_modular<T> || is_runtime_modular<T>;}
constexpr vector<T> multiply_ntt(const vector<T> &a, const vector<T>& b){
    int n = 1;
    while(n < a.size() + b.size()) n *= 2;
    vector<T> _a(n), _b(n);
    copy(begin(a), end(a), begin(_a)); copy(begin(b), end(b), begin(_b));
    FFT::fft(_a); FFT::fft(_b);
    for (int i = 0; i < n; i++) _a[i] *= _b[i];
    FFT::fft(_a, 1);
    return _a;
}
template <Number T> requires requires {requires is_modular<T> || is_runtime_modular<T>;}
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
    if(a.size() + o.size() >= 128){
        if constexpr (S == MulType::ntt){
            return multiply_ntt(a, o);
        }else if constexpr (S == MulType::floating) {
            return multiply_fft_cut(a, o);
        }
    } return multiply_naive(a, o);
}

} //end fft
#undef cpx
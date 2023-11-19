#include <iostream>
#include <vector>
#include <sstream>
using namespace std;

struct Decimal{
    typedef long long ll;
    static constexpr int baselen = 9;
    static constexpr int base = 1e9;
    static constexpr int maxdigit = base - 1;
    alignas(256) static constexpr int pow10[baselen + 1] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};

    vector<int> v;
    int sign;

    explicit operator bool() const{ return sign != 0; }
    static char dig(int a){
        return '0' + a;
    }
    static int put(char a){
        return a - '0';
    }
    operator string() const{
        if(sign == 0) return "0";
        stringstream ans;
        if(sign < 0) ans << '-';
        ans << *v.rbegin();
        for(auto i = v.rbegin()+1; i != v.rend(); ++i){
            string s(baselen, '0');
            for(int j = 0; j < baselen; ++j){
                s[baselen-1-j] = dig(*i % pow10[j + 1] / pow10[j]);
            } ans << s;
        }return ans.str();
    }
    friend ostream& operator<<(ostream& f, const Decimal& o){
        return f << (string)o;
    }
    friend istream& operator>>(istream& f, Decimal& o){
        string s; f >> s; o = s; return f;
    }
    int log10_abs() const {
        if(v.empty()) return 0;
        int s = (v.size() - 1) * baselen + (int)(to_string(v.back()).size());
        return s;
    }
    int len() const {
        if(v.empty()) return 1;
        int s = (v.size() - 1) * baselen;
        s += (sign < 0);
        s += (int)(to_string(v.back()).size());
        return s;
    }
    Decimal(): v(0), sign(0){}
    // <typename T, enable_if_t<is_integral_v<T> && is_signed_v<T>, T*>
    template <signed_integral T>
    Decimal(T n) : v(), sign((n > 0) - (n < 0)){
        if(sign < 0) {
            while (n < 0){
                v.push_back(-(int)(n % base));
                n /= base;
            }
        }else {
            while (n > 0){
                v.push_back((int)(n % base));
                n /= base;
            }
        }
    }
    template <unsigned_integral T>
    Decimal(T n) : v(), sign(n > 0){
        while (n > 0){
            v.push_back((int) (n % base));
            n /= base;
        }
    }
    Decimal(const string &s): v(2 + s.size() / baselen), sign{1}{
        if(s.empty() || s == "0"){v={}, sign=0;}
        else{
            if(s[0] == '-') sign = -1;
            int j = s.size(), d = (sign == -1), i = 0;
            while(j - baselen > d){
                for(int k = j-baselen; k < j; ++k)
                    v[i] += put(s[k]) * pow10[j - 1 - k];
                j -= baselen;
                ++i;
            }for(int k = d; k < j; ++k){v[i] += put(s[k]) * pow10[j - 1 - k];}
            reduce(v);
        }
    }
    static vector<int>& reduce(vector<int>& a){
        while(!a.empty() && !a.back()) {a.pop_back();}
        return a;
    }
    static int get(const vector<int>& v, int idx){
        if(0 <= idx && idx < (int)v.size()) return v[idx];
        return 0;
    };
    static vector<int> add(const vector<int>& a, const vector<int>& b){
        int n = max(a.size(), b.size());
        vector<int> ans(n+1, 0);
        for(int i = 0; i < n; ++i){
            ans[i] += get(a, i) + get(b, i);
            if(ans[i] >= base){
                ans[i] -= base;
                ans[i+1]++;
            }
        } reduce(ans);
        return ans;
    }
    static vector<int> sub(const vector<int>& a, const vector<int>& b, int& sign){
        int n = max(a.size(), b.size());
        int k = n-1;
        while(k >= 0 && get(a, k) == get(b, k)) --k;
        vector<int> ans(k+1);
        if(ans.empty()){sign = 0; return {};}
        if(get(a, k) < get(b, k)){
            sign = -1;
            for(int i = 0; i <= k; ++i){
                ans[i] += get(b, i) - get(a, i);
                if(ans[i] < 0){
                    ans[i] += base;
                    ans[i+1]--;
                }
            }
        }
        else{ sign = 1;
            for(int i = 0; i <= k; ++i){
                ans[i] += get(a, i) - get(b, i);
                if(ans[i] < 0){
                    ans[i] += base;
                    ans[i+1]--;
                }
            }
        }
        reduce(ans);
        return ans;
    }
    static vector<int> mul(const vector<int>& a, const vector<int>& b) {
        if (a.empty() || b.empty()) return {};

        int n = a.size(), m = b.size();
        if (n > 4 && m > 4 && n + m > 32) {
            auto ans = karatsuba_mul(a, b);
            reduce(ans);
            return ans;
        }
        vector<int> ans(n + m + 1);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) {
                long long tmp = (long long) a[i] * b[j];
                ans[i + j] += tmp % base;
                if (ans[i + j] >= base) {
                    ans[i + j] -= base;
                    ans[i + j + 1]++;
                }
                ans[i + j + 1] += tmp / base;
                if (ans[i + j + 1] >= base) {
                    ans[i + j + 1] -= base;
                    ans[i + j + 2]++;
                }
            }
        reduce(ans);
        return ans;
    }
    using ull = unsigned long long;
    static vector<ull> karatsuba_transform(const vector<int>& v){
        int n = ((int)v.size() + 1) / 2 * 3;
        vector<ull> c(n);
        auto at = [&v](int i) -> ull { return (i < (int)v.size()) ? v[i] : 0;};
        for(int i = 0, j = 0; j < n; i += 2, j += 3) {
            c[j] = at(i) % 1000000;
            c[j + 1] = (at(i) / 1000000) + 1000 * (at(i + 1) % 1000);
            c[j + 2] = at(i + 1) / 1000;
        }
        return c;
    }

    static vector<int> karatsuba_reverse(const vector<ull>& v){
        int n = ((int)v.size() + 2) / 3 * 2;
        vector<int> c(n);
        auto at = [&v](int i) -> ull { return (i < (int)v.size()) ? v[i] : 0;};
        for(int i = 0, j = 0; i < n; i += 2, j += 3) {
            c[i] = at(j) + 1000000 * (at(j + 1) % 1000);
            c[i + 1] = (at(j + 1) / 1000) + 1000 * at(j + 2);
        }
        return c;
    }
    static void karatsuba_mul_recursion(const ull* a, const ull* c, ull* inplace, int st){
        if(st <= 3) {
            for (int i = 0; i < (1 << st); ++i)
                for (int j = 0; j < (1 << st); ++j)
                    inplace[i + j] += a[i] * c[j];
            return;
        } const int t = 1 << (st - 1), l = 1 << st;
        karatsuba_mul_recursion(a, c, inplace, st-1);
        karatsuba_mul_recursion(a + t, c + t, inplace + l, st-1);
        vector<ull> rec(l), wrt(l);
        for(int i = 0; i < t; ++i)
            rec[i] = a[i] + a[t + i];
        for(int i = 0; i < t; ++i)
            rec[t + i] = c[i] + c[t + i];
        karatsuba_mul_recursion(rec.data(), rec.data()+t, wrt.data(), st-1);
        for(int i = 0; i < t; ++i){
            ull x = inplace[i+t], y = inplace[i+l];
            inplace[i+t] += wrt[i] - inplace[i] - y;
            inplace[i+l] += wrt[i+t] - inplace[i+t+l] - x;
        }
    }
    static vector<int> karatsuba_mul(const vector<int>& a, const vector<int>& b){
        auto va = karatsuba_transform(a), vb = karatsuba_transform(b);
        int f = 0;
        while((va.size() > (1u << f)) || (vb.size() > (1u << f) )) ++f;
        va.resize(1u << f); vb.resize(1u << f);
        auto vc = vector<ull>(1u << (f+1));
        karatsuba_mul_recursion(va.data(), vb.data(), vc.data(), f);
        for(int i = 0; i + 1 < (1 << (f + 1)); ++i){
            vc[i+1] += vc[i] / 1000000;
            vc[i] %= 1000000;
        }
        return karatsuba_reverse(vc);
    }

    static Decimal addsub(const Decimal& a, const Decimal& b, bool adds = true){
        if (b.v.empty()) return a;
        if (a.v.empty()) {
            Decimal ans = b;
            ans.sign *= (adds ? 1 : -1);
            return ans;
        }
        Decimal ans;
        if ((a.sign == b.sign) == adds){
            ans.v = add(a.v, b.v);
            ans.sign = a.sign;
        } else{
            ans.v = sub(a.v, b.v, ans.sign);
            ans.sign *= a.sign;
        }if(ans.v.empty()) ans.sign = 0;
        return ans;
    }
    Decimal operator+(const Decimal& o) const {
        return addsub(*this, o, true);
    }
    Decimal operator-(const Decimal& o) const {
        return addsub(*this, o, false);
    }
    Decimal operator*(const Decimal& o) const {
        Decimal ans;
        ans.v = mul(v, o.v);
        ans.sign = sign * o.sign;
        return ans;
    }
    auto operator<=>(const Decimal& o) const {
        if(sign != o.sign) return sign <=> o.sign;
        if(sign == 0) return 0 <=> 0;
        if(sign > 0){
            if(v.size() != o.v.size() || v.size() == 0) return v.size() <=> o.v.size();
            for(int i = v.size() - 1; i > 0; --i) if(v[i] != o.v[i]) return v[i] <=> o.v[i];
            return v[0] <=> o.v[0];
        }else {
            if(v.size() != o.v.size() || v.size() == 0) return o.v.size() <=> v.size();
            for(int i = v.size() - 1; i > 0; --i) if(v[i] != o.v[i]) return o.v[i] <=> v[i];
            return o.v[0] <=> v[0];
        }
    }
    int operator%(int mod) const {
        int ans = 0;
        for(int n = v.size(), i = n-1; i >= 0; --i){
            ans = ((long long)ans * Decimal::base + v[i]) % mod;
        } if(sign < 0) ans = -ans;
        if(ans < 0) ans += abs(mod);
        return ans;
    }
    static int idiv(Decimal& a, int b){
        int n = a.v.size(), m = n-1;
        vector<int> ans(n);
        long long now = 0;
        while(m >= 0){
            now = Decimal::base * now + a.v[m];
            ans[m] = now / b;
            now %= b;
            --m;
        }
        for(int i = 0; i + 1 < n; ++i) if(ans[i] >= Decimal::base){
            ans[i] -= Decimal::base;
            ans[i+1]++;
        }
        while(!ans.empty() && ans.back() == 0) ans.pop_back();
        a.v.swap(ans);
        if(a.v.empty()) a.sign = 0;
        return now;
    }
    Decimal operator/(int n){
        Decimal ans = *this;
        idiv(ans, n); return ans;
    }
    /* being lazy.. */
    Decimal& operator+=(const Decimal& o){
        return *this = *this + o;
    }
    Decimal& operator-=(const Decimal& o){
        return *this = *this - o;
    }
    Decimal& operator*=(const Decimal& o){
        return *this = *this + o;
    }
};

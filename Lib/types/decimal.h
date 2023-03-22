#include <iostream>
#include <vector>
using namespace std;

struct Decimal{
    typedef long long ll;
    static constexpr int baselen = 9;
    static constexpr int base = 1e9;
    static constexpr int maxdigit = base - 1;
    alignas(256) static constexpr int pow10[baselen + 1] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};

    vector<int> v;
    int sign;

    operator bool() const{ return sign != 0; }
    operator string() const{
        if(sign == 0) return "0";
        string ans;
        ans.reserve(baselen*v.size() + 1);
        if(sign < 0) ans += '-';
        ans += to_string(*v.rbegin());
        for(auto i = v.rbegin()+1; i != v.rend(); ++i){
            string s(baselen, '0');
            for(int j = 0; j < baselen; ++j){
                s[baselen-1-j] = '0' + (*i % pow10[j+1] / pow10[j]);
            } ans += s;
        }return ans;
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
    // concept:
    // constexpr Decimal(signed_integral n)
    template <typename T, enable_if_t<is_integral_v<T> && is_signed_v<T>, T*> = nullptr>
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
    // concept:
    // constexpr Decimal(unsigned_integral n)
    template <typename T, enable_if_t<is_unsigned_v<T>, T*> = nullptr> 
    Decimal(T n) : v(), sign(n > 0){
        while (n > 0){
            v.push_back((int) (n % base));
            n /= base;
        }
    }
    Decimal(const string &s):v(2+s.size()/baselen), sign{1}{
        if(s.empty() || s == "0"){v={}, sign=0;}
        else{
            if(s[0] == '-') sign = -1;
            int j = s.size(), d = (sign == -1), i = 0;
            while(j - baselen > d){
                for(int k = j-baselen; k < j; ++k){
                    v[i] += (s[k]-'0') * pow10[j-1-k];
                }j -= baselen;
                ++i;
            }for(int k = d; k < j; ++k){v[i] += (s[k]-'0') * pow10[j-1-k];}
            reduce(v);
        }
    }
    static vector<int>& reduce(vector<int>& a){
        while(!a.empty() && !a.back()) a.pop_back();
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
    static vector<int> mul(const vector<int>& a, const vector<int>& b){
        if(a.empty() || b.empty()) return {};
        int n = a.size(), m = b.size();
        vector<int> ans(n+m+1);
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < m; ++j){
                long long tmp = (long long)a[i] * b[j];
                ans[i+j] += tmp % base;
                if(ans[i+j] >= base){
                    ans[i+j] -= base; ans[i+j+1]++;
                }
                ans[i+j+1] += tmp / base;
                if(ans[i+j+1] >= base){
                    ans[i+j+1] -= base; ans[i+j+2]++;
                }
            }
        reduce(ans);
        return ans;
    }
    // static pair<vector<int>, vector<int>> divmod(const vector<int>& a, const vector<int>& b){// todo
    //     assert(!b.empty());  
    // } 

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
        }//if(ans.v.empty()) ans.sign = 0;
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
};

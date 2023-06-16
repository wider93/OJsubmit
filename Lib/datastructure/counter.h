#pragma once
#include "iostream"
#include <vector>
#include <iterator>
using namespace std;
constexpr int maxn = 28;
struct Counter{
    char x[maxn];
    constexpr Counter():x{}{}
    constexpr Counter(int n, int m = 1):x{}{x[n] = m;}
    Counter(const string& s):x{}{
        char y[26]{};
        for(auto i: s) y[i-'A']++;
        for(int i = 0; i < 26; ++i)x[y[i]]++;
        x[0] = 0;
    }
    constexpr Counter operator+(const Counter& o) const {
        Counter ans;
        for(int i = 0; i < maxn; ++i) ans.x[i] = x[i] + o.x[i];
        return ans;
    }
    constexpr Counter& operator+=(const Counter& o) {
        for(int i = 0; i < maxn; ++i) x[i] += o.x[i];
        return *this;
    }
    constexpr Counter operator*(int n) const {
        Counter ans;
        for(int i = 0; i < maxn; ++i) ans.x[i] = x[i] * n;
        return ans;
    }
    constexpr Counter moves(int from, int to) const{
        Counter ans = *this;
        ans.x[from]--; ans.x[from-1]++;
        ans.x[to]--; ans.x[to+1]++;
        ans.x[0] = 0; //in case from = 1
        return ans;
    }
    constexpr auto operator<=> (const Counter& o) const = default;
    friend ostream& operator<<(ostream& f, const Counter& a){
        for(int i = maxn-1; i >= 0; --i) if(a.x[i]) f << i << ": " << (int)a.x[i] <<", ";
        return f << '\n';
    }
};
vector<Counter> dp[maxn][maxn] = {};
vector<Counter> alldp[maxn] = {};
void init(){
    for(int n = 1; n < maxn; ++n) dp[n][1] = {n};
}
vector<Counter> space(int, int);
vector<Counter> allspace(int n, int m){
    if(n == 0) return {{}};
    if(m == 0) return {};
    vector<Counter> ans;
    for(int i = m; i >= 1; --i) {
        auto res = space(n, i);
        ans.insert( ans.end(), make_move_iterator(res.begin()), make_move_iterator(res.end()) );
    } return ans;
}
vector<Counter> allspace(int n){ return allspace(n, n);}
vector<Counter> space(int n, int m){//sum, m is used at least once
    if(m > n) return {};
    if(m == 0 || dp[n][m].size()) return dp[n][m];
    vector<Counter> ans;
    for(int k = n/m; k >= 1; --k){
        auto res = allspace(n-m*k, m - 1);
        for(auto &j: res) j += Counter(m, k);
        ans.insert( ans.end(), make_move_iterator(res.begin()), make_move_iterator(res.end()) );
    }
    return dp[n][m] = ans;
}

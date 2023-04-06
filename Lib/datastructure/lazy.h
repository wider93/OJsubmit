#pragma once
#include "concepts"
#include "vector"
using namespace std;
template<typename T, typename S>
constexpr bool Distributive = true;
template <typename Val, typename Op>
concept AppropriateLazy = requires(const Val& a, const Val& b, const Op& f, const Op& g){
    {Val::e()} -> same_as<Val>;
    {Op::e()} -> same_as<Op>;
    {a + b} -> same_as<Val>;
    {g(a) } -> same_as<Val>;
    {f + g} -> same_as<Op>;
} && Distributive<Val, Op>;

template<typename Val, typename Op>
    requires AppropriateLazy<Val, Op>
struct LazySeg { // 1-indexed
private:
    constexpr static int _lg(int n) {
        int ret = 0; while (n > 1 << ret) ++ret;
        return ret;
    }
    void apply(int i, const Op& f) {
        tree[i] = f(tree[i]);
        if (i < sz) lazy[i] = f + lazy[i];
    }
    void push(int i) {
        apply(i << 1, lazy[i]);
        apply(i << 1 | 1, lazy[i]);
        lazy[i] = Op::e();
    }
    void pull(int i) {
        tree[i] = tree[i << 1] + tree[i << 1 | 1];
    }
public:
    const int n, lg, sz;
    vector<Val> tree;
    vector<Op> lazy;
    constexpr LazySeg(int n) : n(n), lg(_lg(n)), sz(1 << lg), tree(sz << 1, Val::e()), lazy(sz, Op::e()) {}
    constexpr LazySeg(vector<Val>&& inits) : LazySeg(inits.size()){
        copy(begin(inits), end(inits), begin(tree)+sz);
        init();
    }
    void set(int i, const Val& val){tree[i + sz] = val;}
    void init() { for (int i = sz - 1; i > 0; --i) pull(i); }
    void update(int i, const Op& f) {
        i += sz;
        for (int j = lg; j; j--) push(i >> j);
        apply(i, f);
        while(i >>= 1) pull(i);
    }
    void update(int l, int r, const Op& f) {// [l, r)
        l += sz, r += sz;
        for (int i = lg; i; i--) {
            push(l >> i);
            if((l >> i) != (r >> i)) push(r >> i);
        }
        for (int L = l, R = r; L < R; L >>= 1, R >>= 1) {
            if (L & 1) apply(L++, f);
            if (R & 1) apply(--R, f);
        }
        for (int i = 1; i <= lg; i++) {
            push(l >> i);
            if((l >> i) != (r >> i)) push(r >> i);
        }
    }
    Val query(int i) {
        i += sz;
        for (int j = lg; j; j--) push(i >> j);
        return tree[i];
    }
    Val query(int l, int r) {// [l, r)
        l += sz, r += sz ;
        Val L = Val::e(), R = Val::e();
        for (int i = lg; i > 0; i--) {
            push(l >> i);
            if((l >> i) != (r >> i)) push(r >> i);
        }
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) L = L + tree[l++];
            if (r & 1) R = tree[--r] + R;
        }
        return L + R;
    }
};

typedef long long ll;
struct Val {
    ll v; int sz;
    constexpr Val(ll v, int sz = 1):v(v), sz(sz){};
    static constexpr Val e(){return {0, 1};}
    constexpr Val operator+(const Val& o) const{
        return {v + o.v, sz + o.sz};
    }
};

struct Op {
    ll a;
    constexpr Op(ll a): a(a){}
    static constexpr Op e(){ return 0; }
    constexpr Op operator+(const Op& o) const {
        return {a + o.a};
    }
    constexpr Val operator()(const Val& x) const {
        return {x.v + a * x.sz, x.sz};
    }
};

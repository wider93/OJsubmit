#pragma once
#include <vector>
#include <array>
#include <algorithm>
using namespace std;
template <signed_integral T>
struct Box{
    const vector<T> bounds;
    static bool nonempty(const vector<T>& v){return !v.empty() && all_of(v.begin(), v.end(), [](T a){return a != 0;});}
    struct iterator {
        vector<T> f;
        const vector<T>* s;
        iterator(const vector<T>& x, const vector<T>* y): f(x), s(y){}
        auto operator<=>(const iterator& o) const { return f <=> o.f; }
        bool operator==(const iterator& o) const { return f == o.f; }
        iterator& operator++() {
            int n = ssize(*s) - 1;
            while(n > 0 && f[n] == (*s)[n] - 1) --n;
            for(++f[n++]; n < ssize(*s); ++n) f[n] = 0;
            return *this;
        }
        iterator& operator--() {
            int n = ssize(*s) - 1;
            while(n > 0 && f[n] == 0) --n;
            for(--f[n++]; n < ssize(*s); ++n) f[n] = (*s)[n] - 1;
            return *this;
        }
        auto&& operator*() const {
            return this -> f;
        }
    } const b, e;
    static vector<T> make_end(const vector<int>& v, T v0){
        vector<T> ans(v.size());
        if(!ans.empty()) ans[0] = v0;
        return ans;
    }
    explicit Box(const vector<T>& v): bounds{v},
        b{make_end(v, 0), &bounds},
        e{make_end(v, nonempty(v) ? v[0] : T{}), &bounds} {}
    iterator begin() const{ return b; }
    iterator end() const { return e; }
};
template <typename S>
struct Product{
    const vector<vector<S>> elems;
    const vector<int> bounds;
    static bool nonempty(const vector<vector<S>>& v){return all_of(v.begin(), v.end(), [](auto& a){return !a.empty();});}
    struct iterator : Box<int>::iterator{
        const vector<vector<S>> *e;
        iterator(const vector<int>& x, const vector<int>* y, const vector<vector<S>>* z)
            : Box<int>::iterator(x, y), e(z) {}
        vector<S> operator*() const {
            vector<S> substitute(s -> size());
            for(int i = 0; i < ssize(*s); ++i) substitute[i] = (*e)[i][f[i]];
            return substitute;
        }
    } const zero;
    static vector<int> init(const vector<vector<S>>& v){
        vector<int> ans(size(v));
        for(int i = 0; i < ssize(v); ++i) ans[i] = ssize(v[i]);
        return ans;
    }
    explicit Product(const vector<vector<S>> &v) : elems(v), bounds(init(v)), zero{vector(size(v), 0), &bounds, &elems} {}
    iterator begin() const {return zero;}
    iterator end() const {
        auto ans = zero;
        ans.f[0] = (*ans.s)[0];
        return ans;
    }
};

template <signed_integral T, int N>
struct BoxStatic{
    const array<T, N> bounds;
    static bool nonempty(const array<T, N>& v){return all_of(v.begin(), v.end(), [](T a){return a != 0;});}
    struct iterator {
        array<T, N> f;
        const array<T, N>* s;
        constexpr iterator(T v, const array<T, N>* y): f{v}, s(y){}
        constexpr auto operator<=>(const iterator& o) const { return f <=> o.f; }
        constexpr bool operator==(const iterator& o) const { return f == o.f; }
        constexpr iterator& operator++() {
            int n = N - 1;
            while(n > 0 && f[n] == (*s)[n] - 1) --n;
            for(++f[n++]; n < ssize(*s); ++n) f[n] = 0;
            return *this;
        }
        constexpr iterator& operator--() {
            int n = N - 1;
            while(n > 0 && f[n] == 0) --n;
            for(--f[n++]; n < ssize(*s); ++n) f[n] = (*s)[n] - 1;
            return *this;
        }
        constexpr auto&& operator*() const {
            return this -> f;
        }
    } const b, e;
    constexpr explicit BoxStatic(const array<T, N>& v): bounds{v}, b{0, &bounds}, e{nonempty(v) ? v[0] : 0, &bounds} {}
    constexpr iterator begin() const{ return b; }
    constexpr iterator end() const { return e; }
};
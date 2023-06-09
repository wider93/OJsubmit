#include <memory>
#include <vector>

using namespace std;
template <signed_integral T>
struct box{
    const vector<T> bounds;
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
    } const zero;
    explicit box(const vector<T>& v): bounds{v}, zero{vector(ssize(v), 0), &bounds} {}
    iterator begin() const{ return zero; }
    iterator end() const {
        auto ans = zero;
        ans.f[0] = (*ans.s)[0];
        return ans;
    }
};
template <typename S>
struct product{
    const vector<vector<S>> elems;
    const vector<int> bounds;
    struct iterator : box<int>::iterator{
        const vector<vector<S>> *e;
        iterator(const vector<int>& x, const vector<int>* y, const vector<vector<S>>* z)
            : box<int>::iterator(x, y), e(z) {}
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
    explicit product(const vector<vector<S>> &v) : elems(v), bounds(init(v)), zero{vector(ssize(v), 0), &bounds, &elems} {}
    iterator begin() const {return zero;}
    iterator end() const {
        auto ans = zero;
        ans.f[0] = (*ans.s)[0];
        return ans;
    }
};


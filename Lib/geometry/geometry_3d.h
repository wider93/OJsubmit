#include "iostream"
using namespace std;
template <typename F>
struct R3{
    F a, b, c;
    R3(const F& a, const F& b, const F& c): a(a), b(b), c(c){}
    R3(): a{}, b{}, c{}{}
    R3 operator*(const R3& o) const{ return {b*o.c-c*o.b, c*o.a-a*o.c, a*o.b-b*o.a}; }
    R3 operator*(const F &o) const{ return {a * o, b * o, c * o}; }
    R3 operator/(const F &o) const{ return {a / o, b / o, c / o}; }
    R3 operator+(const R3& o) const{ return {a+o.a, b+o.b, c + o.c}; }
    R3 operator-(const R3& o) const{ return {a-o.a, b-o.b, c - o.c}; }
    R3 operator-() const{ return {-a, -b, -c}; }
    F inner(const R3& o) const {return a * o.a + b * o.b + c * o.c;}
    R3& reduce() {
        a.reduce(); b.reduce(); c.reduce();
        return *this;
    }
    friend ostream& operator<<(ostream& f, const R3& o){ return f << o.a << ' ' << o.b << ' ' << o.c << '\n'; }
    friend istream& operator>>(istream& f, R3& o){ return f >> o.a >> o.b >> o.c; }
};
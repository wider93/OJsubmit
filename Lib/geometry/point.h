#pragma once
#include <iostream>
using namespace std;
template <typename T>
struct Point{
    T x, y;
    constexpr Point(T x = 0, T y = 0):x(x), y(y){}
    constexpr Point operator+(const Point& o) const { return {x + o.x, y + o.y}; }
    constexpr Point operator-(const Point& o) const { return {x - o.x, y - o.y}; }
    constexpr Point& operator+=(const Point& o) { x += o.x; y += o.y; return *this; }
    constexpr Point& operator-=(const Point& o) { x -= o.x; y -= o.y; return *this;}
    constexpr auto operator<=>(const Point& o) const = default;
    friend istream& operator>>(istream& f, Point& p){return f >> p.x >> p.y;}
    friend ostream& operator<<(ostream& f, const Point& p){return f << p.x << ' ' << p.y;}
};
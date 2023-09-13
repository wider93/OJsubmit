#pragma once
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
};
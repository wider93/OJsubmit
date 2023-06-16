#pragma once
#include <algorithm>
#include <vector>
#include "types/integers.h"
#include "geometry/point.h"
using namespace std;
template <typename T>
T ccw(const Point<T> &A, const Point<T> &B, const Point<T> &P) {
    return (A.x - P.x) * (B.y - P.y) - (A.y - P.y) * (B.x - P.x);
}
template <integral T>
wider_t<T> ccw(const Point<T> &A, const Point<T> &B, const Point<T> &P) {
    return wider_t<T>(A.x - P.x) * (B.y - P.y) - wider_t<T>(A.y - P.y) * (B.x - P.x);
}
template <typename T>
vector<int> convexHullIndex(const vector<Point<T>> &q){
    int n = q.size();
    if (n <= 2) {
        vector<int> ans(n);
        for(int i = 0; i < n; ++i) ans[i] = i;
        return ans;
    }
    vector<int> p(n);
    for(int i = 0; i < n; ++i) p[i] = i;
    swap(p[0], *min_element(begin(p), end(p), [&q](int a, int b){return q[a] < q[b];}));
    const auto p0 = q[p[0]];
    auto comp = [&p0, &q](int a, int b) {
        auto v = ccw(q[a], q[b], p0);
        if (v) return v > 0;
        return q[a] < q[b];
    };
    sort(p.begin() + 1, p.end(), comp);

    vector<int> s = {p[0], p[1]};
    int now = 2, next = 2;
    while (next < n) {
        while (now >= 2 && ccw(q[s[now-2]], q[s[now-1]], q[p[next]]) <= 0) {
            s.pop_back(); now--;
        }
        s.push_back(p[next++]);
        ++now;
    } return s;
}
template <typename T>
vector<Point<T>> toCoordinates(const vector<Point<T>> &q, const vector<int>& indices){
    vector<Point<T>> ans(indices.size());
    for(int i = 0, n = indices.size(); i < n; ++i)
        ans[i] = q[indices[i]];
    return ans;
}
template <typename T>
vector<Point<T>> convexHull(const vector<Point<T>> &q){
    auto indices = convexHullIndex(q);
    return toCoordinates(q, indices);
}

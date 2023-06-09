#pragma once
#include <iostream>
#include <cassert>
#include <traits.h>
#include "vector"
using namespace std;
void fastio(){
    cin.tie(0) ->sync_with_stdio(0);
}
template <typename T>
concept Printable = requires (T a){
    {cout << a } -> same_as<ostream&>;
};
template <Container T>
void println(T a, ostream& f = cout){
    auto now = a.begin(), end = a.end();
    if(now != end) {
        f << *(now++);
    } while(now != end){
        f << ' ' << *(now++);
    }
}

template <Printable T>
ostream& operator<<(ostream& f, const vector<T> &v){
    println(v, f); return f;
}

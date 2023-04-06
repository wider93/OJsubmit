#pragma once
#include <iostream>
#include <cassert>
#include <traits.h>
using namespace std;

template <Container T>
void println(T a, ostream& f = cout){
    auto now = a.begin(), end = a.end();
    if(now != end) {
        f << *(now++);
    } while(now != end){
        f << ' ' << *(now++);
    }
}
template <Container T>
ostream& operator<<(ostream& f, T v){
    println(v, f); return f;
}
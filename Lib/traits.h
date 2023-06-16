#pragma once
#include "concepts"
using namespace std;

template<typename T>
concept Number = requires(const T& a, T& b, const T& c) {
    {a + c} -> same_as<T>;
    {a - c} -> same_as<T>;
    {a * c} -> same_as<T>;
    {a / c} -> same_as<T>;
    { -a } -> same_as<T>;
    {b += c} -> same_as<T&>;
    {b -= c} -> same_as<T&>;
    {b *= c} -> same_as<T&>;
    {b /= c} -> same_as<T&>;
    { T{ 0 } };	 { T{ 1 } };	
};

template <typename T>
concept Container = requires(const T& a) {
    {a.begin()};
    {a.end()};
};

//template<typename T>
//concept Hashable = requires(T a){
//    { hash<T>{}(a) } -> convertible_to<size_t>;
//};
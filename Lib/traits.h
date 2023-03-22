#pragma once
#include <concepts>
#include <functional>
using namespace std;
template<typename T>
concept Hashable = requires(T a){
    { hash<T>{}(a) } -> convertible_to<size_t>;
};

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

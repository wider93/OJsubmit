#include "iostream"
#include "vector"
#include "algorithm"
#include "cassert"
#include "cmath"
using namespace std;
//template <typename T>
//constexpr array<T, 5> d(T p){
//    T q = T(1) - p;
//    return {p*p*p*p, 4*p*p*p*q, 6*p*p*q*q, 4*p*q*q*q, q*q*q*q};
//}
//#undef NDEBUG
constexpr int l = 11;
auto x(auto i){
    return i+1;
}
int main(){
    cout << fixed;
    long long s = 0;
    cout << x(s) << ' ' << x(1.);
}
#include "iostream"
#include "types/poly.h"
#include "types/modulo.h"
using namespace std;

using Z = Z_n<998244353>;
using poly = Poly<Z, MulType::ntt>;
int main(){
    cin.tie(0) -> sync_with_stdio(0);
    int n; cin >> n;
    vector<int> parent(n+1);
    for(int i = 2; i <= n; ++i) cin >> parent[i];
    vector<Z> p(n+1);
    for(int i = 1; i <= n; ++i) cin >> p[i];
}


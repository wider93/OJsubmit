#include "myfile.cpp"
#include "sub.h"

using namespace std;





int main(){
    cin.tie(0) -> sync_with_stdio(0);
    i32 tc; cin >> tc;
    for(i32 i = 0; i < tc; ++i){
        i64 f; cin >> f;
        auto t = factorize(f);
        cout << t.size() << ' ';
        for(const auto& i: t) cout << i << ' ';
        cout << '\n';
    }
}


#pragma once
#include "vector"
#include "bit"
#include "cassert"
using namespace std;

template<class T>
struct RMQ {
    int n, depth;
	vector<vector<T>> jmp;

	RMQ(const vector<T>& v): n(size(v)), depth(bit_width(unsigned(n + 1))), jmp(depth, v) {
		for(int i = 0; i < depth-1; ++i)
            for(int j = 0; j < n; ++j)
                jmp[i+1][j] = min(jmp[i][j], jmp[i][min(n - 1, j + (1 << i))]);
	};

	T query(int a, int b) const{
		assert(a < b); // or return inf if a == b
		auto dep = bit_width((unsigned)(b - a)) - 1;
		return min(jmp[dep][a], jmp[dep][b - (1 << dep)]);
	}
};
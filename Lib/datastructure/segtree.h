#pragma once
#include "vector"
using namespace std;

template <typename T, T combine(const T&, const T&), T mass_combine(const T&, const T&, int), T identity=0>
struct LazySeg{
	vector<T> tree;
    vector<int> length;
	vector<T> modify;
	int N;
	LazySeg(int n): N(n), tree(2 * n, identity), length(2 * n, 1), modify(2 * n){
        for(int i = n - 1; i >= 1; --i) length[i] = length[2 * i] + length[2 * i + 1];
	}
    LazySeg(int n, const vector<T> &init): LazySeg(n){
        for(int i = n; i < 2 * n; ++i) tree[i] = init[i - n];
        rebuild();
	}
    void pop(int i) {
	    tree[i >> 1] = combine(mass_combine(tree[i],     modify[i >> 1], length[i]    ),
                               mass_combine(tree[i ^ 1], modify[i >> 1], length[i ^ 1]));
	}
    void popUp(int i){
        for (; i > 1; i >>= 1) pop(i);
	}
    void rebuild(){
        for(int i = N - 1; i >= 1; --i) tree[i] = combine(tree[2 * i], tree[2 * i + 1]);
	}
    void singleUpdate(int i, T v){ // update i-th value to v
		tree[i += N] = v;
		while(i >>= 1) tree[i] = combine(tree[i << 1], tree[i << 1 | 1]);
	}
    void modifierHelper(int i,T p){
        tree[i] = mass_combine(tree[i], p, length[i]);
        modify[i] = combine(modify[i], p);
    }
    void rangeUpdate(int l, int r, T p) {
        l += N; r += N;
        int LCopy = l;
        int RCopy = r;
        while(l < r) {
            if (l & 1) { modifierHelper(l++, p); }
            if (r & 1) { modifierHelper(--r, p); }
            l >>= 1, r >>= 1;
        }
        popUp(LCopy);
        popUp(RCopy - 1);
	}

    T realValue(int i) const{
        T v = tree[i];
        for (int j = i >> 1; j > 0; j >>= 1) v = mass_combine(v, modify[j], length[i]);
        return v;
    }

    T rangeQuery(int l, int r) const{ // query for [l, r)
    	l += N, r += N;
    	T rl = identity, rr = identity;
    	while(l < r){
    		if(l & 1) rl = combine(rl, realValue(l++));
    		if(r & 1) rr = combine(realValue(--r), rr);
    		l >>= 1, r >>= 1;
    	} return combine(rl, rr);
    }
};

template <typename T, T combine(const T&, const T&)>
struct Seg{
	const int N;
	const T identity;
	vector<T> tree;
	Seg(int n, const T& iden = T{}):N(n), identity(iden), tree(2 * n, identity){}
	Seg(int n, vector<T> &&init): Seg(n){
		copy(init.begin(), init.end(), tree.begin() + N);
		rebuild();
	}
	void rebuild(){
		for(int i = N - 1; i >= 1; --i) tree[i] = combine(tree[2 * i], tree[2 * i + 1]);
	}
	void singleUpdate(int i, const T& v){ // update i-th value to v
		tree[i += N] = v;
		while(i >>= 1) tree[i] = combine(tree[i << 1], tree[i << 1 | 1]);
	}
	T rangeQuery(int l, int r) const{ // query for [l, r)
		l += N, r += N;
		T rl = identity, rr = identity;
		while(l < r){
			if(l & 1) rl = combine(rl, tree[l++]);
			if(r & 1) rr = combine(tree[(--r)], rr);
			l >>= 1, r >>= 1;}
		return combine(rl, rr);
	}
};
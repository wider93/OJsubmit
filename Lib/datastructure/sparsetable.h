#pragma once
#include "vector"
using namespace std;

struct SparseTable{
    int depth;
    vector<vector<int>> table;
    SparseTable(const vector<int>& f, int depth): depth(depth), table(depth, f){
        for(int i = 1; i < depth; ++i){
            for(int j = 0; j < ssize(f); ++j)
                table[i][j] = table[i - 1][table[i - 1][j]];
        }
    }
    int composite(int nod, int steps) const{
        for(int i = 0; i < depth && steps; ++i, steps>>=1)
            if(steps & 1) nod = table[i][nod];
        return nod;
    }
};
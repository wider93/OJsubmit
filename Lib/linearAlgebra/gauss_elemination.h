#pragma once
#include <vector>
using namespace std;

template <typename T>
void ref(vector<vector<T>>& mat){
    const int m = mat.size(), n = mat[0].size();
    int now = 0;
    for(int i = 0; i < n; ++i){
        bool hasOne = false;
        int o = now;
        for(int k = now; k < m; ++k) if(mat[k][i]) {o = k; hasOne = true; break;}
        if(!hasOne) continue;
        swap(mat[now], mat[o]);
        T const val = mat[now][i];
        for(int k = i; k < n; ++k)
            mat[now][k] /= val;
        for(int j = now + 1; j < m; ++j){
            T const val = mat[j][i];
            for(int k = i; k < n; ++k)
                mat[j][k] -= val * mat[now][k];
        } now++;
    }
}
template <typename T>
void rref(vector<vector<T>>& mat){
    const int m = mat.size(), n = mat[0].size();
    int now = 0;
    for(int i = 0; i < n; ++i){
        bool hasOne = false;
        int o = now;
        for(int k = now; k < m; ++k) if(mat[k][i]) {o = k; hasOne = true; break;}
        if(!hasOne) continue;
        swap(mat[now], mat[o]);
        T const val = mat[now][i];
        for(int k = i; k < n; ++k)
            mat[now][k] /= val;
        for(int j = 0; j < now; ++j){
            T const val = mat[j][i];
            for(int k = i; k < n; ++k)
                mat[j][k] -= val * mat[now][k];
        }
        for(int j = now + 1; j < m; ++j){
            T const val = mat[j][i];
            for(int k = i; k < n; ++k)
                mat[j][k] -= val * mat[now][k];
        } now++;
    }
}
template <typename T>
vector<T> solveEquation(const vector<vector<T>>& mat, int val = 0){ 
    /** used on rref-ed n*(n+k) matrix, and should have unique solutions
    if val == 0, then it calculates for first variable (a.k.a n-th  col)
    else it calculates for (n+val) -th col
    **/
    vector<T> ans(mat.size());
    val += mat.size();
    for(int i = mat.size() - 1; i >= 0; --i){
        T x = mat[i][val];
        for(int j = mat.size() - 1; j > i; --j){
            x -= mat[i][j] * ans[j];
        } ans[i] = x;
    } return ans;
}

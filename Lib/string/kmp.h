#pragma once
#include "vector"
#include "ranges"
using namespace std;

template<random_access_iterator It>
struct kmp{
    int n;
    using T = remove_reference_t<decltype (* ((It)nullptr))>;
    vector<T> needle;
    vector<int> fail;
    kmp(It b, It e):n(e - b), needle(b, e), fail(n){
        for(int j = 0, i = 1; i < n; ++i){
            while(j > 0 && b[i] != b[j]) j = fail[j-1];
            if (b[i] == b[j]) ++j, fail[i] = j;
            else fail[i] = 0;
        }
    }
    template<random_access_iterator It2>
    requires (same_as<T, remove_reference_t<decltype (* ((It2)nullptr))>>)
    vector<int> search(It2 b, It2 e){
        int m = distance(b, e) - n + 1;
        vector<int> ans;
		for(int j = 0; b < e; ++b){
			while (j > 0 && !(j < n && needle[j] == *b))
				j = fail[j-1];
			if(j < n && needle[j] == *b){
				if (j + 1 == n) ans.push_back(m - distance(b, e));
			}++j;
		} return ans;
    }
};
#pragma once
#include <tuple>
#include "datastructure/RMQ.h"
#include "algorithm"
#include "datastructure/sparsetable.h"
struct LCA {
	vector<int> time;
	vector<long long> dist;
	RMQ<pair<int, int>> rmq;

	explicit LCA(const vector<vector<pair<int, int>>>& graph) : time(size(graph), -(1 << 30)), dist(size(graph)), rmq(dfs(graph)) {}

	vector<pair<int, int>> dfs(const vector<vector<pair<int, int>>>& graph) {
		vector<tuple<int, int, int, long long>> q(1);
		vector<pair<int, int>> ret;
		int T = 0;
		while (!q.empty()) {
			auto [v, p, d, di] = q.back();
			q.pop_back();
			if (d) ret.emplace_back(d, p);
			time[v] = T++;
			dist[v] = di;
			for(const auto &e: graph[v]) if (e.first != p)
				q.emplace_back(e.first, v, d+1, di + e.second);
		}
		return ret;
	}
	int query(int a, int b) const {
		if (a == b) return a;
		a = time[a], b = time[b];
		return rmq.query(min(a, b), max(a, b)).second;
	}
	long long distance(int a, int b) const{
		int lca = query(a, b);
		return dist[a] + dist[b] - 2 * dist[lca];
	}
};

namespace graph{
    template<typename T> int getNode(const T& a);
    template<> int getNode(const int& a){return a;}
    template<> int getNode(const pair<int, int>& a){return a.first;}
}
struct LCAjump : SparseTable {
    vector<int> fromRoot;
    template<class T>
    static constexpr pair<vector<int>, vector<int>> forestToParent(const vector<vector<T>>& g, vector<int> roots = {0}){
        vector<int> ans(g.size(), -1);
        vector<int> dist(g.size());
        for(auto &i: roots) ans[i] = i;
        while(!roots.empty()){
            vector<int> rep;
            for(auto i: roots){
                for(auto &j: g[i]){
                    int jn = graph::getNode(j);
                    if(ans[jn] != -1) continue;
                    ans[jn] = i;
                    dist[jn] = dist[i] + 1;
                    rep.push_back(jn);
                }
            } roots.swap(rep);
        } return {ans, dist};
    }
    LCAjump(const vector<int> &f, vector<int>&& dist) :
            SparseTable(f, bit_width((unsigned)*max_element(dist.begin(), dist.end())) + 1u), fromRoot(dist) {}

    template<class NodeType>
    LCAjump(const vector<vector<NodeType>> &graph, vector<int> roots = {0}) :
        LCAjump(forestToParent(graph, roots).first,forestToParent(graph, roots).second){ }
    int query(int a, int b) const{
        if (fromRoot[a] < fromRoot[b]) swap(a, b);
        a = composite(a, fromRoot[a] - fromRoot[b]);
        if (a == b) return a;
        for (int i = depth-1; i--;) {
            int c = table[i][a], d = table[i][b];
            if (c != d) a = c, b = d;
        } return table[0][a];
    }
    int distance(int a, int b) const{
		int lca = query(a, b);
		return fromRoot[a] + fromRoot[b] - 2 * fromRoot[lca];
	}
};
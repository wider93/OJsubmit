#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <queue>
// #include "demo.cpp"

using namespace std;

int main(){
    auto isBracketHeader = [&](const string& s)-> bool{
        return s.size() > 2 and s[0] == '<' and s.back() == '>';
    };
    auto isHeader = [&](const string& s)-> bool{
        return s.size() > 2 and ((s[0] == '<' and s.back() == '>') or (s[0] == '"' and s.back() == '"'));
    };
    auto headerToFilename = [&](const string &s){
        return s.substr(1, s.size() - 2);
    };
    
    auto stl = [&](const string& fileName) -> unordered_set<string>{
        auto target = ifstream(fileName);
        stringstream f;
        f << target.rdbuf();
        unordered_set<string> ans;
        string s;
        while(f >> s){
            if(isBracketHeader(s)) 
                ans.insert(headerToFilename(s));
        } return ans;
    }("c++headers.txt");
    
    // for(auto &s : stl) cout << "#include " << s << '\n';
        cout << '\n';
    auto getHeaders = [&](string fileName) -> vector<string>{
        if(isHeader(fileName)) fileName = headerToFilename(fileName);
        if(stl.contains(fileName)) return {};
        vector<string> ans;
        auto f = ifstream(fileName);
        string s;
        while(f >> s){
            if(s.size() < 2) continue;
            else if(s.substr(0, 2) == "/*") {
                while(s.size() >= 2 && s.substr(s.size() - 2, 2) != "*/" && f >> s);
                continue;
            }
            else if(s.substr(0, 2) == "//") {
                getline(f, s);
                continue;
            }
            if(s == "#include") {
                f >> s;
                assert(isHeader(s));
                ans.push_back(s);
            } else if(s.size() > 10 and s.substr(0, 8) == "#include"){
                auto t = s.substr(8, s.size() - 8);
                assert(isHeader(t));
                ans.push_back(t);
            } else {
                ;
            }
        } return ans;
    };
    for(auto i: getHeaders("demo.cpp")) cout << i << '\n';

    
    
    auto alignHeaders = [&](const string& source) -> string{
        // for [a, b] in ans: a placed ahead of b
        queue<string> q;
        unordered_set<string> r; r.reserve(80);
        unordered_map<string, vector<string>> prereq; prereq.reserve(80);
        unordered_map<string, int> left; left.reserve(80);
        q.push(source);
        r.insert(source);
        while(!q.empty()){
            auto s = q.front();
            q.pop();
            if(stl.contains(s)) continue;
            for(const auto &t: getHeaders(s)){
                // cout << s << ", "<< t << '\n';
                prereq[s].push_back(t);
                left[t]++;
                if(r.contains(t)) continue;
                r.insert(t);
                q.push(t);
            }
        }
        vector<string> head; head.reserve(r.size());
        q.push(source); // note source not pushed in head
        vector<string> files; files.reserve(r.size());
        while(!q.empty()){
            const auto &t = q.front();
            for(const auto &i: prereq[t]){
                left[i]--;
                if(!left[i]){
                    q.push(i);
                    if(isBracketHeader(i)) head.push_back(i);
                    else files.push_back(i);
                }
            }q.pop();
        }
        stringstream ans;
        for (auto &i: head) ans << "#include " << i << '\n';
        ans << "using namespace std;\n";
        for (auto &i: files) ans << "contents of {" << i << "}\n";
        ans << "contents of {" << source << "}\n";
        return ans.str();
    };
    cout << '\n';

    auto merged = alignHeaders("demo.cpp");
    cout << merged;
}

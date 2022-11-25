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


bool isBracketHeader(const string& s){
    return s.size() > 2 and s[0] == '<' and s.back() == '>';
};
bool isHeader(const string& s){
    return s.size() > 2 and ((s[0] == '<' and s.back() == '>') or (s[0] == '"' and s.back() == '"'));
};
string headerToFilename(const string &s){
    return s.substr(1, s.size() - 2);
};
bool isBlank(const string& s){
    for(auto &i: s) if(i != ' ' && i != '\n' && i != '\t') return false;
    return true;
}

auto stl = [](const string& fileName) -> unordered_set<string>{
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
bool usingNamespaceStd(const string &s){
    return s.find("using") != string::npos && s.find("namespace") != string::npos && s.find("std;") != string::npos;
}
struct CppFile{
    vector<string> headers;
    string contents;
    CppFile():headers{}, contents{}{}
    CppFile(string fileName):headers{}, contents{}{
        if(isHeader(fileName)) fileName = headerToFilename(fileName);
        if(stl.contains(fileName)) return;
        auto f = ifstream(fileName);
        string s;
        stringstream contentsStream;
        while(getline(f, s) && !usingNamespaceStd(s)){
            if(s.size() > 10 and s.substr(0, 8) == "#include"){
                size_t i = 8; while(i < s.size() && s[i] == ' ') ++i;
                auto t = s.substr(i, s.size() - i);
                assert(isHeader(t));
                headers.push_back(t);
            } else if (!isBlank(s)){
                contentsStream << s << '\n';
            }
        }
        while(getline(f, s)){
            if(isBlank(s)) continue;
            contentsStream << s << '\n';
        }contents = contentsStream.str();
    }
};
namespace std{
    template <> struct hash<CppFile>{
        size_t operator()(CppFile const& s) const noexcept{
            return hash{}(s.contents) ^ s.headers.size();
        }
    };
}
int main(int argv, char* argc[]){
    auto makeFile = [&](const string& source) -> string{
        // for [a, b] in ans: a placed ahead of b
        queue<string> q;
        unordered_set<string> r; r.reserve(818);
        unordered_map<string, CppFile> prereq; prereq.reserve(818);
        unordered_map<string, int> left; left.reserve(818);
        q.push(source);
        r.insert(source);
        while(!q.empty()){
            auto s = q.front();
            q.pop();
            if(stl.contains(s)) continue;
            prereq[s] = CppFile(s);
            auto &file = prereq[s];
            cout << s << endl;
            for(const auto &t: file.headers){
                left[t]++;
                if(r.contains(t)) continue;
                r.insert(t);
                q.push(t);
            }
        }
        vector<string> head; head.reserve(r.size());
        q.push(source); 
        vector<string> files; files.reserve(r.size());
        while(!q.empty()){
            const auto &t = q.front();
            if(isBracketHeader(t)) head.push_back(t);
            else files.push_back(prereq[t].contents);
            for(const auto &i: prereq[t].headers){
                left[i]--;
                if(!left[i]){
                    q.push(i);
                }
            }q.pop();
        }
        stringstream ans;
        for (auto &i: head) 
            ans << "#include " << i << '\n';
        ans << "\nusing namespace std;\n";
        for (auto i = files.rbegin(); i != files.rend(); ++i) 
            ans << *i << '\n';
        return ans.str();
    };

    string s = "demo.cpp";
    if(argv != 1) s = argc[1];
    auto merged = makeFile(s);
    cout << merged;
}


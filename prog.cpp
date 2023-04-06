#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-interfaces-global-init"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ranges>
#include <vector>
#include <sstream>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>

using namespace std;

constexpr char ONCE[] = "pragma once";

bool isHeader(const string& s){
    return s.size() > 2 and ((s[0] == '<' and s.back() == '>') or (s[0] == '"' and s.back() == '"'));
};
string headerToFilename(const string &s){
    return s.substr(1, s.size() - 2);
};
constexpr bool isBlank(char i){
    return i == ' ' || i == '\n' || i == '\t' || i == '\r';
};
constexpr bool isAllBlank(const string& s){
    return ranges::all_of(s.begin(), s.end(), isBlank);
}
bool isBlankOrComment(const string& s){
    for(int j = 0, n = s.size(); j < n; ++j){
        if(isBlank(s[j])) continue;
        return j + 1 < n && s[j] == '/' && s[j+1] == '/';
    } return true;
}

const unordered_set<string> stl = [](const string& fileName){
    auto target = ifstream(fileName);
    stringstream f;
    f << target.rdbuf();
    unordered_set<string> ans;
    string s;
    while(f >> s){
        if(isHeader(s))
            ans.insert(headerToFilename(s));
    } return ans;
}("../c++headers.txt");
bool isSTL(const string& s){
    return stl.count(s);
};
bool usingNamespaceStd(const string &s){
    auto a = s.find("using"), b = s.find("namespace"), c = s.find("std;");
    return a < b && b < c && c < string::npos;
}
struct CppFile{
    static constexpr string_view prefixes[]{"./", "../", "../Lib/"};
    vector<string> headers;
    string contents;
    CppFile():headers{}, contents{}{}
    CppFile(string fileName):headers{}, contents{}{
        if(isHeader(fileName)) fileName = headerToFilename(fileName);
        if(stl.contains(fileName)) return;
        fileName = firstOccurence(fileName);
        auto f = ifstream(fileName);
        string s;
        stringstream contentsStream;
        while(getline(f, s) && !usingNamespaceStd(s)){
            if(s.size() > 10 and s.substr(0, 8) == "#include"){
                size_t i = 8; while(i < s.size() && s[i] == ' ') ++i;
                auto t = s.substr(i, s.size() - i);
                assert(isHeader(t));
                headers.push_back(t.substr(1, t.size() - 2));
            }else if (!(isBlankOrComment(s) || s.find(ONCE) != string::npos)){
                contentsStream << s << '\n';
            }
        }
        while(getline(f, s)){
            if(isBlankOrComment(s) || s.find(ONCE) != string::npos) continue;
            contentsStream << s << '\n';
        }contents = contentsStream.str();
    }
    static string firstOccurence(string& fileName){
        for(auto &prefix: prefixes){
            auto trial = string(prefix.begin(), prefix.end()) + fileName;
            auto f = ifstream(trial);
            string s;
            if(f >> s) {
                return trial;
            }
        } return "";
    }
};
namespace std{
    template <> struct hash<CppFile>{
        size_t operator()(CppFile const& s) const noexcept{
            return hash{}(s.contents) ^ s.headers.size();
        }
    };
}

string makeFile(const string& source){
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
        if(isSTL(t)) head.push_back(t);
        else files.push_back(prereq[t].contents);
        for(const auto &i: prereq[t].headers){
            left[i]--;
            if(!left[i])
                q.push(i);
        }q.pop();
    }
    stringstream ans;
    for (auto &i: head)
        ans << "#include <" << i << ">\n";
    ans << "\nusing namespace std;\n";
    for (auto &file : std::ranges::reverse_view(files)) {
        if (file.empty()) continue;
        ans << file;
    } return ans.str();
};
int main(int argv, char* argc[]){
    string s = "main.cpp";
//    for(auto &i: stl) cout << i << ' ';
//    cout << '\mod';
    if(argv != 1) s = argc[1];
    auto merged = makeFile(s);
    cout << merged;
}


#pragma clang diagnostic pop
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
// #include "demo.cpp"

using namespace std;

int main(){
    auto stl = [](const string& filename) -> vector<string>{
        auto target = ifstream(filename);
        stringstream f;
        f << target.rdbuf();
        vector<string> ans;
        string s;
        while(f >> s){
            if(s[0] == '<' && s.back() == '>') ans.push_back(s);
        } return ans;
    }("c++headers.txt");
    for(auto &s : stl) cout << "#include " << s << '\n';
}
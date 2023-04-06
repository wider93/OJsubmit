#include <fstream>
#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argv, char* argc[]){
    string s = argv > 1 ? argc[1] : "main.cpp";
    //auto t = R"(dir && .\prog.exe )"+ s +" > ../source.cpp"s;
    auto t = R"(.\prog.exe )"+ s +" > ../source.cpp"s;
    cout << t << '\n';
    system(t.c_str());
}
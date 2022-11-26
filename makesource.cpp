#include <fstream>
#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argv, char* argc[]){
    string s = argv > 1 ? argc[1] : "demo.cpp";
    auto t = R"(cd .. && .\cmake-build-debug\prog.exe )"+ s +" > source.cpp"s;
    cout << t << '\n';
    system(t.c_str());
}
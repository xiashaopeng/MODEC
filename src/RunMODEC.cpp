//#include <iostream>
#include "ModecClass.h"

using namespace std;

ModecClass MODEC;
int main(int argc, char *argv[]) {
    //std::ios::sync_with_stdio(false);
    //cin.tie(NULL);
    setlocale(LC_ALL, "chs");
    MODEC.RunModec(argc, argv);
    InfoMessage::ends.push_back(clock());
    InfoMessage::EndInfo();
    return 0;
}

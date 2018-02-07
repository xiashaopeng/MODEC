#ifndef _ERROR_MESSAGE_H
#define _ERROR_MESSAGE_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#ifdef _WIN32||defined _WIN64
#include <direct.h>
#else
#include <limits.h>
#include <unistd.h>
#endif
#include <time.h>
using namespace std;

class InfoMessage {
  public:
    static clock_t start;
    static vector<clock_t> ends;
    static string InputName;
    static string error_info;
    static vector<string> warning_info;
    static void ErrorMessage(string message, int mode);
    static void EndInfo();
};


#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <assert.h>
#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <cctype>
#include <cstdio>
#include <stdio.h>
#include <float.h>

using namespace std;
//0 for Hairpin
//1 for internal
//2 for bulige
vector<vector<double>> energy;

vector<vector<double>> W;
vector<vector<double>> V;
vector<vector<double>> M;
vector<vector<double>> M1;
vector<vector<int>> trace_W;
vector<vector<int>> trace_M1;
vector<vector<int>> trace_V;
int TracebackV(int i,int j);
int FindhlforTracebackV(int i, int j);
vector<char> set_array;
int Trace;
string model;
int a;
int b;
double c;

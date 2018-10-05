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
#include <math.h>
#include <stdio.h>

using namespace std;
vector<vector<double>> TracebackM(int i ,int j,string model1,string model2);
vector<vector<double>> TracebackX(int i ,int j,string model1,string model2);
vector<vector<double>> TracebackY(int i ,int j,string model1,string model2);
vector<vector<double>> TracebackH(int i ,int j,string model1,string model2);
vector<vector<double>> TracebackP(int i ,int j,string model1,string model2);
vector<vector<double>> TracebackQ(int i ,int j,string model1,string model2);
vector<vector<double>> TraceForNeedlemanWunschGotoh(string model1,string model2);
void TestForNeedlemanWunschGotoh(string model1,string model2);
void NeedlemanWunschGotoh(string model1,string model2,int i, int j);
void NeedlemanWunschGotohBackward(string model1,string model2);
void SmithWatermanGotoh(string model1,string model2);
using namespace std;
vector<vector<double>> Trace;
vector<vector<double>> sum;
vector<vector<double>> M;
vector<vector<double>> X;
vector<vector<double>> Y;
vector<vector<double>> bM;
vector<vector<double>> bX;
vector<vector<double>> bY;
vector<vector<double>>H;
vector<vector<double>>P;
vector<vector<double>>Q;
int d;
int e;

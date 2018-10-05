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
#include <float.h> // FLT_EPSILON, DBL_EPSILON, LDBL_EPSILON
#include <math.h>

using namespace std;
vector<vector<float>> TracebackVm(int i ,int j,string model1,string model2);
vector<vector<float>> TracebackVx(int i ,int j,string model1,string model2);
vector<vector<float>> TracebackVy(int i ,int j,string model1,string model2);
vector<vector<float>> TracebackH(int i ,int j,string model1,string model2);
vector<vector<float>> TracebackP(int i ,int j,string model1,string model2);
vector<vector<float>> TracebackQ(int i ,int j,string model1,string model2);
vector<vector<float>> TraceForNeedlemanWunschGotoh(string model1,string model2);
float MaxOfVm(int i,int j);
using namespace std;
vector<vector<float>> Trace;
vector<vector<float>> Fm;
vector<vector<float>> Fx;
vector<vector<float>> Fy;
vector<vector<float>> Vm;
vector<vector<float>> Vx;
vector<vector<float>> Vy;
vector<vector<float>>Bm;
vector<vector<float>>Bx;
vector<vector<float>>By;
vector<vector<int>>trace_Vm;//Vm.0 Vx.1 Vy.2
vector<vector<float>> xiyi_percentage;
vector<char> array1;
vector<char> array2;
float alpha;
float beta;
float delta;
float qxi;
float qyj;
float xy_percentage_;
float xypi_percentage_;

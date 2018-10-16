#include "main.h"

int FormBasePair(char a,char b){
  if(a=='G'||b=='G'){
    if(a=='U'||b=='U'||a=='C'||b=='C')
      return 1;
    return 0;
  }
  if(a=='A'||b=='A'){
    if(a=='U'||b=='U')
      return 1;
    return 0;
  }
  return 0;
}

double StackingScore(int i, int j){
  if(model[i+1]=='U'||model[j-1]=='U'){
    if(model[i]=='U'||model[j]=='U')
      return -0.5;
    return -2.0;
  }
  if(model[i]=='U'||model[j]=='U')
    return -2.0;
  return -3.0;
}


double F1ForHairpin(int i,int j){
  return energy[0][j-i-1];
}

double F2(int i,int j,int h,int l){
  if(h==i+1&&l==j-1)
    return StackingScore(i, j);
  if(l==j-1)
    return energy[2][h-i-1];
  if(h==i+1)
    return energy[2][j-l-1];
  return energy[1][h-i-1+j-l-1];
}

double MinOfF2PlusVForVij(int i, int j){
  double min_score=INFINITY;
  for(int h=i+1;h<j-1;h++){
    for(int l=h+1;l<j;l++){
      if(F2(i,j,h,l)+V[h][l]<min_score){
	min_score=F2(i,j,h,l)+V[h][l];
      }
    }
  }
  return min_score;
}


double FindVij(int i,int j){
  int trace_V_candidate1;
  double min_candidate1;
  if(!FormBasePair(model[i],model[j]))
    return INFINITY;
  if(F1ForHairpin(i,j)>MinOfF2PlusVForVij(i,j)){
    trace_V_candidate1=2;
    min_candidate1=MinOfF2PlusVForVij(i,j);
  }
  else{
    trace_V_candidate1=1;
    min_candidate1=F1ForHairpin(i,j);
  }
  if(min_candidate1>(M[i+1][j-1]+a+b)){
    trace_V[i][j]=3;
    return M[i+1][j-1]+a+b;
  }
  trace_V[i][j]=trace_V_candidate1;
  return min_candidate1;
}

//For M
double MinOfM1PlusM1(int i, int j){
  double min_score=INFINITY;
  for(int k=i;k<j;k++){
    if(M1[i][k]+M1[k+1][j]<min_score){
      min_score=M1[i][k]+M1[k+1][j];
    }
  }
  return min_score;
}

double FindKFromM1PlusM1ForM(int i, int j){
  double min_score=INFINITY;
  int min_k;
  for(int k=i;k<j;k++){
    if(M1[i][k]+M1[k+1][j]<min_score){
      min_score=M1[i][k]+M1[k+1][j];
      min_k=k;
    }
  }
  return min_k;
}

//for M1
double FindM1ij(int i,int j){
  int trace_M1_candidate1;
  int trace_M1_candidate2;
  double min_candidate1;
  double min_candidate2;
  
  if(M[i][j]>(V[i][j]+b)){
    trace_M1_candidate1=2;
    min_candidate1=V[i][j]+b;
  }
  else{
    trace_M1_candidate1=1;
    min_candidate1=M[i][j];
  }
  if((M1[i+1][j]+c)>(M1[i][j-1]+c)){
    trace_M1_candidate2=4;
    min_candidate2=M1[i][j-1]+c;
  }
  else{
    trace_M1_candidate2=3;
    min_candidate2=M1[i+1][j]+c;
  }
  if(min_candidate1>min_candidate2){
    trace_M1[i][j]=trace_M1_candidate2;
    return min_candidate2;
  }
  trace_M1[i][j]=trace_M1_candidate1;
  return min_candidate1;
}

//For W
double MinOfWPlusW(int i, int j){
  double min_score=INFINITY;
  for(int k=i;k<j;k++){
    if(W[i][k]+W[k+1][j]<min_score){
      min_score=W[i][k]+W[k+1][j];
    }
  }
  return min_score;
}
double FindKFromWPlusW(int i, int j){
  double min_score=INFINITY;
  int min_k=0;
  for(int k=i;k<j;k++){
    if(W[i][k]+W[k+1][j]<min_score){
      min_score=W[i][k]+W[k+1][j];
      min_k=k;
    }
  }
  return min_k;
}
double FindWij(int i,int j){
  int trace_W_candidate1;
  int trace_W_candidate2;
  double min_candidate1;
  double min_candidate2;
  
  if(W[i+1][j]>W[i][j-1]){
    trace_W_candidate1=2;
    min_candidate1=W[i][j-1];
  }
  else{
    trace_W_candidate1=1;
    min_candidate1=W[i+1][j];
  }
  if(V[i][j]>MinOfWPlusW(i,j)){
    trace_W_candidate2=4;
    min_candidate2=MinOfWPlusW(i,j);
  }
  else{
    trace_W_candidate2=3;
    min_candidate2=V[i][j];
  }
  if(min_candidate1>min_candidate2){
    trace_W[i][j]=trace_W_candidate2;
    return min_candidate2;
  }
  trace_W[i][j]=trace_W_candidate1;
  return min_candidate1;
}

void Zukker(string model){
  trace_W.resize(model.size()+1,vector<int>(model.size()+1));
  trace_V.resize(model.size()+1,vector<int>(model.size()+1));
  trace_M1.resize(model.size()+1,vector<int>(model.size()+1));
  W.resize(model.size()+1,vector<double>(model.size()+1));
  V.resize(model.size()+1,vector<double>(model.size()+1));
  M.resize(model.size()+1,vector<double>(model.size()+1));
  M1.resize(model.size()+1,vector<double>(model.size()+1));
  
  for (int i=0;i<model.size();i++) {
    W[i][i]=0;
    V[i][i]=INFINITY;
    M[i][i]=INFINITY;
    M1[i][i]=INFINITY;
  }
  for (int i=1;i<model.size();i++) {
    W[i][i-1]=0;
    V[i][i-1]=INFINITY;
    M[i][i-1]=INFINITY;
    M1[i][i-1]=INFINITY;
  }
  for (int i=0;i<model.size();i++) {
    V[i][i+1]=INFINITY;
    M[i][i+1]=MinOfM1PlusM1(i,i+1);
    W[i][i+1]=FindWij(i,i+1);
    M1[i][i+1]=FindM1ij(i,i+1);
  }
  for(int k=2;k<model.size();k++){
    for(int i=0;i<model.size()-k;i++){
      int j=i+k;
      V[i][j]=FindVij(i,j);
      M[i][j]=MinOfM1PlusM1(i,j);
      W[i][j]=FindWij(i,j);
      M1[i][j]=FindM1ij(i,j);
    }
  }
  cout<<"min_score"<<W[0][model.size()-1]<<endl;
}

void InsertInitialValue(string file){
  ifstream ifs(file);
  if(!ifs.fail()){
    string char_number;
    ifs>>char_number;
    energy.clear();
    for (int i=0;i<3;i++) {
      vector<double> numbers_line_;
      for (int j=0;j<31;j++) {
	if(char_number=="Inf"){
	  numbers_line_.push_back(INFINITY);
	}
	else{
	  numbers_line_.push_back(stod(char_number));
	}
	ifs>>char_number;
      }
      energy.push_back(numbers_line_);
    }
  }
}
	

void TracebackInner(string model){
  vector<vector<int>> trace_array;int i;int j;
  set_array.resize(model.size(),'-');
  int L=model.size()-1;
  trace_array.push_back({0,L});
  while(trace_array.size()){
    i=trace_array.back()[0];
    j=trace_array.back()[1];
    trace_array.pop_back();
    if(i>j-1)
      continue;
    if(trace_W[i][j]==1){
      trace_array.push_back({i+1,j});
      continue;
    }
    if(trace_W[i][j]==2){
     trace_array.push_back({i,j-1});
    continue;
    }
    if(trace_W[i][j]==3){
     int h= TracebackV(i,j);
     continue;
    }
    if(trace_W[i][j]==4){
    int k=FindKFromWPlusW(i, j);
    trace_array.push_back({k+1,j});
    trace_array.push_back({i,k});
    continue;
    }
  }
  for(int j=0;j<model.size();j++){
  cout<<set_array[j];
  }
  cout<<endl;
}

int TracebackM1(int i,int j){
  if(i>j-1)
    return Trace;
  if(trace_M1[i][j]==1){
    int k=FindKFromWPlusW(i, j);
    TracebackM1(i+1,k);
    TracebackM1(k+1,j-1);
  }
  if(trace_M1[i][j]==2){
    TracebackV(i,j);
  }
  if(trace_M1[i][j]==3){
    TracebackM1(i+1,j);
  }
  if(trace_M1[i][j]==4){
    TracebackM1(i,j-1);
  }
}

int TracebackV(int i, int j){
  set_array[i]='(';
  set_array[j]=')';
  if(i>j-1)
    return Trace;
  if(trace_V[i][j]==1){
    return Trace;
  }
  if(trace_V[i][j]==2){
    FindhlforTracebackV(i, j);
  }
  if(trace_V[i][j]==3){
    int k=FindKFromM1PlusM1ForM(i+1,j-1);
    TracebackM1(i+1,k);
    TracebackM1(k+1,j-1);
  }
}
int FindhlforTracebackV(int i, int j){
  double min_score=INFINITY;
  int h_min;
  int l_min;
  for(int h=i+1;h<j-1;h++){
    for(int l=h+1;l<j;l++){
      if(F2(i,j,h,l)+V[h][l]<min_score){
	min_score=F2(i,j,h,l)+V[h][l];
	h_min=h;
	l_min=l;
      }
    }
  }
    return TracebackV(h_min,l_min);
}

//For test
void RunForTest(){
  model="GGGGUUU";
   Zukker(model);
   //test for V-1,2
   assert(V[0][6]==4.1);//F(i0,j6,h1,l5)+V(1,5)
   //test for M1-234
   assert(M1[0][4]<3.51&&M1[0][4]>3.499);//M(1,4)+0.1(c)=3.5
   
}
void RunForKadai(){
  model="CCACCUCGCCGCGAAGAUUGGAUGUAAUGUUCACUUAAGCCUAGGCGUUCGCAAAGAAAUACCUAUGGUAAUUUGACUAGCGGUAACAAUGAAAGAAAAGGUGUAGGCGAGGUGGG";
    Zukker(model);
    TracebackInner(model);
}

int main(void){
    a =	6.0;
    b =	-1.0;
    c= 	0.1;
    InsertInitialValue("sample.fasta");
    RunForKadai();
    //RunForTest();
}



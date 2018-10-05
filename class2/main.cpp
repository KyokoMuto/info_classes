
#include "main.h"

double FindScore(char a,char b){
	return a==b?1:-1;
}

double MaxOfMBackwardMatrix(int i,int j,char a,char b){
	if(bM[i+1][j+1]>bX[i+1][j+1]){
		return bM[i+1][j+1]>bY[i+1][j+1]?bM[i+1][j+1]:bY[i+1][j+1];
	}
	return bX[i+1][j+1]>bY[i+1][j+1]?bX[i+1][j+1]:bY[i+1][j+1];
}
double MaxOfH(int i,int j,char a,char b){
	double maxof2_1=H[i-1][j-1]+FindScore(a,b)>P[i][j]?H[i-1][j-1]+FindScore(a,b):P[i][j];
	double maxof2_2=Q[i][j]>0?Q[i][j]:0;
	return maxof2_1>maxof2_2?maxof2_1:maxof2_2;
}

//For Needleman forward
double MaxOfMMatrix(int i,int j){
	if(M[i-1][j-1]>X[i-1][j-1]){
		return M[i-1][j-1]>Y[i-1][j-1]?M[i-1][j-1]:Y[i-1][j-1];
	}
	return X[i-1][j-1]>Y[i-1][j-1]?X[i-1][j-1]:Y[i-1][j-1];
}
void NeedlemanWunschGotoh(string model1,string model2,int number1,int number2){
	int d=number1;int e=number2;
	//Set initial for M
	for (int i=0;i<model2.size()+1;i++) {
		vector<double> M_line_;
		M_line_.push_back(-1000);
		M.push_back(M_line_);
	}
	M[0].resize(model1.size()+1);
	fill(M[0].begin(),M[0].end(),-1000);
	M[0][0]=0;
	//Set initial for X
	vector<double> X_line_;
	X_line_.push_back(-1000);
	X.push_back(X_line_);
	for (int y=1;y<model1.size()+1;y++) {
		X[0].push_back(-d-(y-1)*e);
	}
	for (int y=1;y<model2.size()+1;y++) {
		vector<double> X_line_;
		X_line_.push_back(-1000);
		X.push_back(X_line_);
	}
	//Set initial for Y
	vector<double> Y_line0_(model1.size()+1);
	fill(Y_line0_.begin(),Y_line0_.end(),-1000);
	Y.push_back(Y_line0_);
	for(int y=1;y<model2.size()+1;y++){
		vector<double> Y_line_;
		Y_line_.push_back(-d-(y-1)*e);
		Y.push_back(Y_line_);
	}
	for(int i=1;i< model2.size()+1;i++){
		for(int j=1;j<model1.size()+1;j++){
			X[i].push_back((M[i-1][j]-d)>(X[i-1][j]-e)?(M[i-1][j]-d):(X[i-1][j]-e));
			Y[i].push_back((M[i][j-1]-d)>(Y[i][j-1]-e)?(M[i][j-1]-d):(Y[i][j-1]-e));
			M[i].push_back(MaxOfMMatrix(i,j)+FindScore(model1[j-1],model2[i-1]));
		}

	}
}
vector<vector<double>> TracebackM(int i,int j,string model1, string model2){
	cout<<model1[j-1]<<" "<<model2[i-1]<<endl;
	if(i==0||j==0)
		return Trace;
	if(M[i][j]==M[i-1][j-1]+FindScore(model1[j-1],model2[i-1]))
		return TracebackM(i-1, j-1,model1,model2);
	if(M[i][j]==X[i-1][j-1]+FindScore(model1[j-1],model2[i-1]))
		return TracebackX(i-1, j-1,model1,model2);
	if(M[i][j]==Y[i-1][j-1]+FindScore(model1[j-1],model2[i-1]))
		return TracebackY(i-1, j-1,model1,model2);
}

vector<vector<double>>TracebackX(int i,int j,string model1, string model2){
	cout<<"-"<<" "<<model2[i-1]<<endl;
	if(i==0||j==0){
		if(j==0){
			for (;i>0;i--){
				cout<<"-"<<" "<<model2[i-1]<<endl;
			}
		}
		for (;j>0;j--){
			cout<<model1[j-1]<<" "<<"-"<<endl;
		}
		return Trace;
	}
	if(X[i][j]==M[i-1][j]-d)
		return TracebackM(i-1,j,model1,model2);
	if(X[i][j]==X[i-1][j]-e)
		return TracebackX(i-1,j,model1,model2);
}
vector<vector<double>> TracebackY(int i, int j,string model1, string model2){
	if(i==0||j==0){
		if(j==0){
			for (;i>0;i--){
				cout<<"-"<<" "<<model2[i-1]<<endl;
			}
		}
		for (;j>0;j--){
			cout<<model1[j-1]<<" "<<"-"<<endl;
		}
		return Trace;
	}
	cout<<model1[j-1]<<" "<<"-"<<endl;
	if(Y[i][j]==M[i][j-1]-d)
		return TracebackM(i,j-1,model1,model2);
	if(Y[i][j]==Y[i][j-1]-e)
		return TracebackY(i,j-1,model1,model2);
}

vector<vector<double>> TraceForNeedlemanWunschGotoh(string model1,string model2,int number1,int number2){
	d=number1;e=number2;
	int j=model1.size();int i=model2.size();
	if(M[i][j]>X[i][j]){
		return M[i][j]>Y[i][j]?TracebackM(i,j,model1,model2):TracebackY(i,j,model1,model2);
	}
	return X[i][j]>Y[i][j]?TracebackX(i,j,model1,model2):TracebackY(i,j,model1,model2);
}

//for smith trace
void SmithWatermanGotoh(string model1,string model2){
	int d=5;int e=1;
	P.resize(model2.size()+1,vector<double>(model1.size()+1));
	Q.resize(model2.size()+1,vector<double>(model1.size()+1));
	H.resize(model2.size()+1,vector<double>(model1.size()+1));
	for (int i=0;i<model2.size()+1;i++) {
		H[i][0]=0;
		P[i][0]=-1000;
		Q[i][0]=-1000;
	}
	for (int j=0;j<model1.size()+1;j++) {
		H[0][j]=0;
		P[0][j]=-1000;
		Q[0][j]=-1000;
	}
	//TODO Delete it
	int u=e;
	int v=d;
	for(int i=1;i<model2.size()+1;i++){
		for(int j=1;j<model1.size()+1;j++){

			P[i][j]=(H[i-1][j]-v)>P[i-1][j]?(H[i-1][j]-v):(P[i-1][j]-u);
			Q[i][j]=(H[i][j-1]-v)>Q[i][j-1]?(H[i][j-1]-v):(Q[i][j-1]-u);
			H[i][j]=MaxOfH(i,j,model1[j-1],model2[i-1]);
		}
	}
}

vector<vector<double>> TracebackH(int i,int j,string model1,string model2){
	if(H[i][j]==0||i==0||j==0){
		return Trace;
	}
	if(H[i][j]==H[i-1][j-1]+FindScore(model1[j-1],model2[i-1])){
		cout<<model1[j-1]<<" "<<model2[i-1]<<endl;
		return TracebackH(i-1, j-1,model1,model2);
	}
	if(H[i][j]==P[i][j])
		return TracebackP(i, j,model1,model2);
	if(H[i][j]==Q[i][j])
		return TracebackQ(i, j,model1,model2);
}

vector<vector<double>>TracebackP(int i,int j,string model1, string model2){
	cout<<"-"<<" "<<model2[i-1]<<endl;
	if(i==0||j==0)
		return Trace;
	if(P[i][j]==H[i-1][j]-d)
		return TracebackH(i-1,j,model1,model2);
	if(P[i][j]==P[i-1][j]-e)
		return TracebackP(i-1,j,model1,model2);
}

vector<vector<double>> TracebackQ(int i, int j,string model1, string model2){
	cout<<model1[j-1]<<" "<<"-"<<endl;
	if(i==0||j==0)
		return Trace;
	if(Q[i][j]==H[i][j-1]-d)
		return TracebackH(i,j-1,model1,model2);
	if(Q[i][j]==Q[i][j-1]-e)
		return TracebackQ(i,j-1,model1,model2);
}

vector<vector<double>> TraceForSmith(string model1,string model2){
	d=5;e=1;
	int maxVal = -1000;
	int max_x=0;
	int max_y=0;
	for(int j = 0; j <model2.size(); j++) {
		for(int i = 0; i <model1.size(); i++) {
			if(H[i][j] > maxVal){
				maxVal = H[i][j];
				max_x=i;
				max_y=j;
			}
		}
	}
	return TracebackH(max_x,max_y,model1,model2);
}

void NeedlemanWunschGotohBackward(string model1,string model2){
	d=5;e=1;
	bM.resize(model2.size()+1,vector<double>(model1.size()+1,0));
	bX.resize(model2.size()+1,vector<double>(model1.size()+1,0));
	bY.resize(model2.size()+1,vector<double>(model1.size()+1,0));
	bM[model2.size()][model1.size()]=0;
	bX[model2.size()][model1.size()]=0;
	bY[model2.size()][model1.size()]=0;
	int k=e;
	int m1_size=model1.size();
	int m2_size=model2.size();
	for (int i=0;i<model2.size();i++) {
		bM[i][model1.size()]=-1000;
		bY[i][model1.size()]=-1000;
		bX[i][model1.size()]=-d-(m2_size-i-2)*k;//(model2.size()-i)*e;
	}
	for (int j=0;j<model1.size();j++) {
		bM[model2.size()][j]=-1000;
		bX[model2.size()][j]=-1000;
		bY[model2.size()][j]=-d-(m1_size-j-2)*k;//(model1.size()-j)*e;
	}
	for(int i=model2.size()-1;i>-1;i--){
		for(int j=model1.size()-1;j>-1;j--){

			bX[i][j]=(bM[i+1][j]-d)>(bX[i+1][j]-e)?(bM[i+1][j]-d):(bX[i+1][j]-e);
			bY[i][j]=(bM[i][j+1]-d)>(bY[i][j+1]-e)?(bM[i][j+1]-d):(bY[i][j+1]-e);
			bM[i][j]=MaxOfMBackwardMatrix(i,j,model1[j],model2[i])+FindScore(model1[j],model2[i]);
		}
	}
}

//For test
void TestForNeedlemanWunschGotoh(string model1,string model2){
	NeedlemanWunschGotoh(model1,model2,2,1);
	assert(X[0][1]==-2&&Y[1][0]==-2);
	assert(X[0][2]==-3&&Y[2][0]==-3);
	assert(M[5][5]==-1&&M[5][6]==1);
	assert(Y[5][6]==-3&&X[5][6]==-4);
}
void RunForKadai1(){
	//TestForNeedlemanWunschGotoh("ACCAGT","ACAGC");
	//vector<vector<double>> test=TraceForNeedlemanWunschGotoh("ACCAGT","ACAGC",2,1);
	string model1="GGAGTGAGGGGAGCAGTTGGGCTGAAGATGGTCAACGCCGAGGGAACGGTAAAGGCGACGGAGCTGTGGCAGACCTGGCTTCCTAACCACGTCCCGTGTTTTGCGGCTCCGCGAGGACTG";
	string model2="CGCATGCGGAGTGAGGGGAGCAGTTGGGAACAGATGGTCCCGCCGAGGGACCGGTGGGCAGACGGGGCCAGCTGTGGCAGACACTGGCTTCTAACCACCGAACGTTCTTTCCGCTCCGC";
	NeedlemanWunschGotoh(model1,model2,5,1);
	vector<vector<double>> f=TraceForNeedlemanWunschGotoh(model1,model2,5,1);
}
void RunForKadai2(){
	string model1="GGAGTGAGGGGAGCAGTTGGGCTGAAGATGGTCAACGCCGAGGGAACGGTAAAGGCGACGGAGCTGTGGCAGACCTGGCTTCCTAACCACGTCCCGTGTTTTGCGGCTCCGCGAGGACTG";
	string model2="CGCATGCGGAGTGAGGGGAGCAGTTGGGAACAGATGGTCCCGCCGAGGGACCGGTGGGCAGACGGGGCCAGCTGTGGCAGACACTGGCTTCTAACCACCGAACGTTCTTTCCGCTCCGC";
	SmithWatermanGotoh(model1,model2);
	TraceForSmith(model1,model2);
}
void RunForKadai3(){
	string model1="GGAGTGAGGGGAGCAGTTGGGCTGAAGATGGTCAACGCCGAGGGAACGGTAAAGGCGACGGAGCTGTGGCAGACCTGGCTTCCTAACCACGTCCCGTGTTTTGCGGCTCCGCGAGGACTG";
	string model2="CGCATGCGGAGTGAGGGGAGCAGTTGGGAACAGATGGTCCCGCCGAGGGACCGGTGGGCAGACGGGGCCAGCTGTGGCAGACACTGGCTTCTAACCACCGAACGTTCTTTCCGCTCCGC";

	NeedlemanWunschGotoh(model1,model2,5,1);
	NeedlemanWunschGotohBackward(model1,model2);
	sum.resize(model2.size()+1,vector<double>(model1.size()+1,0));
	d=5;e=1;
	for(int i=0;i<model2.size()+1;i++){
		for(int j=0;j<model1.size()+1;j++){
			sum[i][j]=M[i][j]+bM[i][j];
			cout<<sum[i][j]<<" ";
		}
		cout<<endl;
	}
}

int main(int kadai_number){
	if(kadai_number==1)
		RunForKadai1();
	if(kadai_number==2)
		RunForKadai2();
	if(kadai_number==3)
		RunForKadai3();
}

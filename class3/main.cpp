#include "main.h"
float PercentageofMatch(char a,char b){
	return a==b?0.125:(0.125/3);
}

//For PairHMM Viterbi
float MaxOfVm(int i,int j){
	if(((1-delta-2*alpha)*Vm[i-1][j-1])>((1-beta-delta)*Vx[i-1][j-1])){
		if(((1-2*alpha-delta)*Vm[i-1][j-1])>((1-beta-delta)*Vy[i-1][j-1])){
			trace_Vm[i][j]=0;
			return (1-2*alpha-delta)*Vm[i-1][j-1];
		}
		trace_Vm[i][j]=2;
		return((1-beta-delta)*Vy[i-1][j-1]);
	}
	if(((1-beta-delta)*Vx[i-1][j-1])>((1-beta-delta)*Vy[i-1][j-1])){
		trace_Vm[i][j]=1;
		return ((1-beta-delta)*Vx[i-1][j-1]);
	}
	trace_Vm[i][j]=2;
	return ((1-beta-delta)*Vy[i-1][j-1]);
}
void PairHMMViterbi(string model1,string model2){
	trace_Vm.resize(model2.size()+1,vector<int>(model1.size()+1));
	Vm.resize(model2.size()+1,vector<float>(model1.size()+1));
	Vx.resize(model2.size()+1,vector<float>(model1.size()+1));
	Vy.resize(model2.size()+1,vector<float>(model1.size()+1));
	Vm[0][0]=1;Vm[1][0]=0;Vm[0][1]=0;Vm[1][1]=0;
	Vx[0][0]=0;Vx[1][0]=alpha*qxi;
	Vy[0][0]=0;Vy[0][1]=alpha*qyj;
	//ToDO takecare of zero
	for (int i=2;i<model2.size()+1;i++) {
		Vm[i][0]=0;
		Vx[i][0]=beta*qxi*Vx[i-1][0];
		Vy[i][0]=0;
	}
	for (int j=2;j<model1.size()+1;j++) {
		Vm[0][j]=0;
		Vy[0][j]=beta*qxi*Vy[0][j-1];
		Vx[0][j]=0;
	}
	for(int i=1;i<model2.size()+1;i++){
		for(int j=1;j<model1.size()+1;j++){
			Vx[i][j]=(alpha*Vm[i-1][j])>(beta*Vx[i-1][j])?(qxi*alpha*Vm[i-1][j]):(qxi*beta*Vx[i-1][j]);
			Vy[i][j]=(alpha*Vm[i][j-1])>(beta*Vy[i][j-1])?(qyj*alpha*Vm[i][j-1]):(qyj*beta*Vy[i][j-1]);
			Vm[i][j]=PercentageofMatch(model1[j-1],model2[i-1])*MaxOfVm(i,j);
		}
	}
}

vector<vector<float>> TracebackVm(int i,int j,string model1, string model2){
	if(i==0||j==0){
		if(i==0){
			for(;j>0;j--){
				array1.push_back(model1[j-1]);
				array2.push_back('-');
			}
		}
		else if(j==0){
			for(;i>0;i--){
				array1.push_back('-');
				array2.push_back(model2[i-1]);
			}
		}
		return Trace;
	}
	array1.push_back(model1[j-1]);
	array2.push_back(model2[i-1]);
	if(trace_Vm[i][j]==0)
		return TracebackVm(i-1, j-1,model1,model2);
	if(Vm[i][j]==PercentageofMatch(model1[j],model2[i])*(1-beta-delta)*Vx[i-1][j-1])
		return TracebackVx(i-1, j-1,model1,model2);
	if(Vm[i][j]==PercentageofMatch(model1[j],model2[i])*(1-beta-delta)*Vy[i-1][j-1])
		return TracebackVy(i-1, j-1,model1,model2);
}

vector<vector<float>>TracebackVx(int i,int j,string model1, string model2){
	if(i==0||j==0){
		if(i==0){
			for(;j>0;j--){
				array1.push_back(model1[j-1]);
				array2.push_back('-');
			}
		}
		else if(j==0){
			for(;i>0;i--){
				array1.push_back('-');
				array2.push_back(model2[i-1]);
			}
		}
		return Trace;
	}
	char arr1='-';
	array1.push_back(arr1);
	array2.push_back(model2[i-1]);
	if(Vx[i][j]==qxi*alpha*Vm[i-1][j])
		return TracebackVm(i-1,j,model1,model2);
	if(Vx[i][j]==qxi*beta*Vx[i-1][j])
		return TracebackVx(i-1,j,model1,model2);
}

vector<vector<float>> TracebackVy(int i, int j,string model1, string model2){
	if(i==0||j==0){
		if(i==0){
			for(;j>0;j--){
				array1.push_back(model1[j-1]);
				array2.push_back('-');
			}
		}
		else if(j==0){
			for(;i>0;i--){
				array1.push_back('-');
				array2.push_back(model2[i-1]);
			}
		}
		return Trace;
	}
	array1.push_back(model1[j-1]);
	array2.push_back('-');
	if(Vy[i][j]==qyj*alpha*Vm[i][j-1])
		return TracebackVm(i,j-1,model1,model2);
	if(Vy[i][j]==qyj*beta*Vy[i][j-1])
		return TracebackVy(i,j-1,model1,model2);
}

vector<vector<float>> TraceForPairHMMViterbi(string model1,string model2){
	int j=model1.size();int i=model2.size();
	if(Vm[i][j]>Vx[i][j]){
		if(Vm[i][j]>Vy[i][j]){
			xypi_percentage_=Vm[i][j];
			return TracebackVm(i,j,model1,model2);
		}
		xypi_percentage_=Vy[i][j];
		return TracebackVy(i,j,model1,model2);
	}
	if(Vx[i][j]>Vy[i][j]){
		xypi_percentage_=Vx[i][j];
		return TracebackVx(i,j,model1,model2);
	}
	xypi_percentage_=Vy[i][j];
	return TracebackVy(i,j,model1,model2);
}

//Pair forward HMM
void PairHMMForward(string model1,string model2){
	Fm.resize(model2.size()+1,vector<float>(model1.size()+1,0));
	Fx.resize(model2.size()+1,vector<float>(model1.size()+1,0));
	Fy.resize(model2.size()+1,vector<float>(model1.size()+1,0));
	Fm[0][0]=1;Fm[1][0]=0;Fm[0][1]=0;Fm[1][1]=0;
	Fx[0][0]=0;Fx[1][0]=alpha*qxi;
	Fy[0][0]=0;Fy[0][1]=alpha*qyj;
	for (int i=2;i<model2.size()+1;i++) {
		Fm[i][0]=0;
		Fx[i][0]=beta*qxi*Fx[i-1][0];
		Fy[i][0]=0;
	}
	for (int j=2;j<model1.size()+1;j++) {
		Fm[0][j]=0;
		Fy[0][j]=beta*qxi*Fy[0][j-1];
		Fx[0][j]=0;
	}

	for(int i=1;i<model2.size()+1;i++){
		for(int j=1;j<model1.size()+1;j++){
			Fm[i][j]=PercentageofMatch(model1[j-1],model2[i-1])*((1-2*alpha-delta)*Fm[i-1][j-1]+(1-beta-delta)*Fx[i-1][j-1]+(1-beta-delta)*Fy[i-1][j-1]);
			Fx[i][j]=qxi*(alpha*Fm[i-1][j]+beta*Fx[i-1][j]);
			Fy[i][j]=qyj*(alpha*Fm[i][j-1]+beta*Fy[i][j-1]);
		}
	}

	xy_percentage_=Fm[model2.size()][model1.size()]+Fx[model2.size()][model1.size()]+Fy[model2.size()][model1.size()];
	cout<<"FE="<<xy_percentage_*delta<<endl;
	cout<<"percentage"<<xypi_percentage_/xy_percentage_<<endl;
}


void PairHMMBackward(string model1,string model2){
	Bm.resize(model2.size()+1,vector<float>(model1.size()+1));
	Bx.resize(model2.size()+1,vector<float>(model1.size()+1));
	By.resize(model2.size()+1,vector<float>(model1.size()+1));
	Bm[model2.size()][model1.size()]=delta;
	Bx[model2.size()][model1.size()]=delta;
	By[model2.size()][model1.size()]=delta;

	for (int i=model2.size()-1;i>0;i--) {
		Bx[i][model1.size()]=beta*qxi*Bx[i+1][model1.size()];
		Bm[i][model1.size()]=alpha*qxi*Bx[i+1][model1.size()];
		By[i][model1.size()]=0;
	}
	for (int j=model1.size()-1;j>0;j--) {
		By[model2.size()][j]=beta*qyj*By[model2.size()][j+1];
		Bm[model2.size()][j]=alpha*qyj*By[model2.size()][j+1];
		Bx[model2.size()][j]=0;
	}
	for(int i=model2.size()-1;i>0;i--){
		for(int j=model1.size()-1;j>0;j--){
			Bm[i][j]=PercentageofMatch(model1[j],model2[i])*(1-2*alpha-delta)*Bm[i+1][j+1]+alpha*qxi*Bx[i+1][j]+alpha*qyj*By[i][j+1];
			Bx[i][j]=PercentageofMatch(model1[j],model2[i])*(1-beta-delta)*Bm[i+1][j+1]+beta*qxi*Bx[i+1][j];
			By[i][j]=PercentageofMatch(model1[j],model2[i])*(1-beta-delta)*Bm[i+1][j+1]+beta*qyj*By[i][j+1];
		}
	}
	xiyi_percentage.resize(model2.size()+1,vector<float>(model1.size()+1,0));
	for(int i=1;i<model2.size()+1;i++){
		for(int j=1;j<model1.size()+1;j++){
			xiyi_percentage[i][j]=(Fm[i][j]*Bm[i][j])/(xy_percentage_*delta);
		}
	}
	for(int i=1;i<model2.size()+1;i++){
		for(int j=1;j<model1.size()+1;j++){
			cout<<" "<<xiyi_percentage[i][j];
		}
		cout<<endl;
	}
}

//For test
void RunForKadai1(){
	alpha=0.1;
	beta=0.2;
	delta=0.1/3;
	qxi=0.25;
	qyj=0.25;
	string model1="CCAGAGCTGTGGCAGACAGTGGCT";
	string model2="CCAGCTGTGCAGACACTGGCTT";
	PairHMMViterbi(model1,model2);
	TraceForPairHMMViterbi(model1,model2);
	for(int j=array1.size()-1;j>-1;j--){
		cout<<array1[j];
	}
	cout<<endl;
	for(int i=array2.size()-1;i>-1;i--){
		cout<<array2[i];
	}
	cout<<endl;
}
void RunForKadai2(){
	alpha=0.1;
	beta=0.2;
	delta=0.1/3;
	qxi=0.25;
	qyj=0.25;
	string model1="CCAGAGCTGTGGCAGACAGTGGCT";ZZ
	string model2="CCAGCTGTGCAGACACTGGCTT";
	PairHMMForward(model1,model2);
}
void RunForKadai3(){
	alpha=0.1;
	beta=0.2;
	delta=0.1/3;
	qxi=0.25;
	qyj=0.25;
	string model1="CCAGAGCTGTGGCAGACAGTGGCT";
	string model2="CCAGCTGTGCAGACACTGGCTT";

	PairHMMForward(model1,model2);
	PairHMMBackward(model1,model2);
}

int main(void){
	RunForKadai1();
	RunForKadai2();
	RunForKadai3();
}



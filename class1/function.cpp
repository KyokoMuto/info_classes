#include "function.h"
int alpha_size_;
vector<char> alpha_;
int state_size;
vector<double> forward_scale_;
vector<vector<int>> hidenn_state_candidates_;
vector<vector<double>> hidenn_state_percentages_;
vector<vector<double>> state_transition_matrix;
vector<vector<double>> output_probability_matrix;
vector<vector<double>> forward_matrix_;
vector<vector<double>> backward_matrix_;
vector<vector<double>> estimation_move_;
vector<vector<double>> estimation_output_;
vector<string> sequence_;

double FindOutPutPercentage(char alpha,int state){
	for (int i=0;i<alpha_size_;i++) {
		if(alpha==alpha_[i]){
			return output_probability_matrix.at(state).at(i);
		}
	}
	return -1;
}

void SetInitialForViterbi(string model){
	hidenn_state_candidates_.clear();
	hidenn_state_percentages_.clear();
	vector<double> first_state_transition_matrix=*state_transition_matrix.begin();
	for (int i=0;i<state_size-1;i++) {
		vector<double> hidenn_state_percentages_line_;
		vector<int> hidenn_state_candidates_line_;
		hidenn_state_candidates_line_.push_back(0);
		double init =log(FindOutPutPercentage(model[0],i))+
			log(first_state_transition_matrix.at(i+1));
		hidenn_state_percentages_line_.push_back(init);
		hidenn_state_percentages_.push_back(hidenn_state_percentages_line_);
		hidenn_state_candidates_.push_back(hidenn_state_candidates_line_);
	}
}

void CalculateAllStatesInViterbi(string model){
	for (int model_number=1;model_number<model.size();model_number++){
		for (int after=1;after<state_size;after++) {
			int hidenn_state_max_candidate_=0;
			double hidenn_state_max_percentage_=-10000;
			for(int before=1;before<state_size;before++) {
				double state_percentage=
					hidenn_state_percentages_.at(before-1).at(model_number-1)
					+log(state_transition_matrix.at(before).at(after))
					+log(FindOutPutPercentage(model[model_number],after-1));
				if(hidenn_state_max_percentage_<state_percentage){
					hidenn_state_max_percentage_=state_percentage;
					hidenn_state_max_candidate_=before;
				}
			}
			hidenn_state_percentages_.at(after-1).push_back(hidenn_state_max_percentage_);
			hidenn_state_candidates_.at(after-1).push_back(hidenn_state_max_candidate_);
		}
	}
}

void OutputHiddenStateForViterbi(string model){
	vector<int> hidden_state_;
	double hidden_state_max_percentage_=-10000;
	int hidden_state_max_candidate_=0;
	for(int i=1;i<state_size;i++){
		double hidden_state_max_percentage_candidate_
			=hidenn_state_percentages_.at(i-1).back();
		if(hidden_state_max_percentage_<hidden_state_max_percentage_candidate_){
			hidden_state_max_percentage_=hidden_state_max_percentage_candidate_;
			hidden_state_max_candidate_=i;
		}
	}
	int before_state=hidden_state_max_candidate_;
	hidden_state_.push_back(before_state);
	before_state=hidenn_state_candidates_.at(before_state-1).back();
	hidden_state_.push_back(before_state);
	for (int model_number=model.size()-2;-1 <model_number;model_number--){
		before_state=hidenn_state_candidates_.at(before_state-1).at(model_number);
		hidden_state_.push_back(before_state);
	}
	while(hidden_state_.size()){
		cout<<hidden_state_.back();
		hidden_state_.pop_back();
	}
}

void InsertInitialValue(string file){
	ifstream ifs(file);
	if(!ifs.fail()){
		char alpha;
		alpha_.clear();
		ifs>>alpha_size_;
		for(int i=0;i<alpha_size_;i++){
			ifs>>alpha;
			alpha_.push_back(alpha);
		}
		ifs>>state_size;
		double state_transition_;
		double output_probability_;
		state_transition_matrix.clear();
		output_probability_matrix.clear();
		for (int i=0;i<state_size;i++) {
			vector<double> state_transition_matrix_line_;
			for (int j=0;j<state_size;j++) {
				ifs>>state_transition_;
				state_transition_matrix_line_.push_back(state_transition_);
			}
			state_transition_matrix.push_back(state_transition_matrix_line_);
		}
		for (int i=0;i<state_size-1;i++) {
			vector<double> output_probability_matrix_line_;
			for (int j=0;j<alpha_size_;j++) {
				ifs>>output_probability_;
				output_probability_matrix_line_.push_back(output_probability_);
			}
			output_probability_matrix.push_back(output_probability_matrix_line_);
		}
	}
}

void SetInitialForForward(string model){
	forward_matrix_.clear();
	forward_scale_.resize(model.size());
	fill(forward_scale_.begin(), forward_scale_.end(), 0);
	vector<double> first_state_transition_matrix=*state_transition_matrix.begin();
	for (int i=0;i<state_size-1;i++) {
		vector<double> forward_line_;
		double element=first_state_transition_matrix.at(i+1)*FindOutPutPercentage(model[0],i);
		forward_scale_[0]+=element;
		forward_line_.push_back(element);
		forward_matrix_.push_back(forward_line_);
	}
	for (int i=0;i<state_size-1;i++) {
		forward_matrix_.at(i).at(0)*=1/forward_scale_[0];
	}
}

void CalculateForForward(string model){
	for (int model_number=1;model_number<model.size();model_number++){
		for (int after=1;after<state_size;after++) {
			double calculated_forward_=0;
			for(int before=1;before<state_size;before++) {
				double element=forward_matrix_.at(before-1).at(model_number-1)*state_transition_matrix.at(before).at(after)*FindOutPutPercentage(model[model_number],after-1);
				calculated_forward_+=element;
			}
			forward_matrix_.at(after-1).push_back(calculated_forward_);
			forward_scale_[model_number]+=calculated_forward_;
		}
		for (int i=0;i<state_size-1;i++) {
			forward_matrix_.at(i).at(model_number)*=1/forward_scale_[model_number];
		}	
	}
}

void SetInitialForBackward(string model){
	backward_matrix_.clear();
	for (int before=1;before<state_size;before++) {
		vector<double> backward_line_;
		double calculated_backward_=0;
		for(int after=1;after<state_size;after++) {
			double element=state_transition_matrix.at(before).at(after)*FindOutPutPercentage(model[model.size()-1],after-1);
			calculated_backward_+=element;
		}
		backward_line_.push_back(calculated_backward_);
		backward_matrix_.push_back(backward_line_);
	}	
	for(int i= 0;i<state_size-1;i++){
		backward_matrix_.at(i).at(0)*=1/forward_scale_[model.size()-1];
	}
}

void CalculateForBackward(string model){
	for (int model_number=model.size()-2;0<model_number;model_number--){
		for (int before=1;before<state_size;before++) {
			double calculated_backward_=0;
			for(int after=1;after<state_size;after++) {
				double element=backward_matrix_.at(after-1).at(model.size()-model_number-2)*state_transition_matrix.at(before).at(after)*FindOutPutPercentage(model[model_number],after-1);
				calculated_backward_+=element;
			}
			backward_matrix_.at(before-1).push_back(calculated_backward_);
		}

		for(int i= 0;i<state_size-1;i++){
			backward_matrix_.at(i).at(model.size()-model_number-1)*=1/forward_scale_[model_number];
		}
	}
}

double SumForward(string model){
	double forward_all_state_=0;
	for(int i=1;i<state_size;i++){
		forward_all_state_+=forward_matrix_.at(i-1).back();
	}
	return forward_all_state_;
}

int FindNumber(char alpha){
	for (int i=0;i<alpha_size_;i++) {
		if(alpha==alpha_[i])
			return i;
	}
	return -1;
}

void Recalculate(string model){
	estimation_move_.resize(state_size, vector<double>(state_size));
	for(int i=0;i<state_size;i++){
		fill(estimation_move_[i].begin(), estimation_move_[i].end(), 0);
	}
	estimation_output_.resize(state_size, vector<double>(alpha_size_));
	for(int i=0;i<state_size;i++){
		fill(estimation_output_[i].begin(), estimation_output_[i].end(), 0);
	}
	//EStep(Estimation) ignore zero...
	for (int model_number=0;model_number<model.size()-2;model_number++){
		for (int before=1;before<state_size;before++) {
			for(int after=1;after<state_size;after++) {
				double element=forward_matrix_.at(before-1).at(model_number)*backward_matrix_.at(after-1).at(model.size()-model_number-3)*state_transition_matrix.at(before).at(after)*FindOutPutPercentage(model[model_number+1],after-1);
				estimation_move_[before][after]+=element;

			}
		}
	}
	for (int model_number=0;model_number<model.size()-1;model_number++){
		for (int before=1;before<state_size;before++) {
			double element1=forward_matrix_.at(before-1).at(model_number)*backward_matrix_.at(before-1).at(model.size()-model_number-2);
			estimation_output_[before][FindNumber(model[model_number])]+=element1;
		}
	}


	//MStep
	vector<double> estimation_output_sum_(state_size,0);
	for(int i=1;i<state_size;i++){
		for(int j=0;j<alpha_size_;j++){
			estimation_output_sum_[i]+=estimation_output_[i][j];
		}
		for(int j=0;j<alpha_size_;j++){
			//Check the percentage of after from alpha and restore outputrate
			output_probability_matrix[i-1][j]=estimation_output_[i][j]/estimation_output_sum_[i];
		}
	}
	vector<double> estimation_move_sum_(state_size,0);
	for(int i=1;i<state_size;i++){
		for(int j=1;j<state_size;j++){
			estimation_move_sum_[i]+=estimation_move_[i][j];
		}
	}
	for(int i=1;i<state_size;i++){
		for(int j=1;j<state_size;j++){
			//Check the percentage of after from before and restore outputrate
			state_transition_matrix[i][j]=estimation_move_[i][j]/estimation_move_sum_[i];
		}
	}
}

//Not used function
void GetSeqences(string file){
	string str1;string str2;
	ifstream ifs1("input.fasta");
	if(!ifs1.fail()){
		getline(ifs1, str1);
		while(ifs1){
			if(str1[0] != '>') {
				if(str2.size()){
					string str3=str1+str2;
					sequence_.push_back(str3);
					getline(ifs1, str1);
					str2="";
					continue;
				}
				getline(ifs1, str2);
			}
			else{
				getline(ifs1, str1);
			}
		}
	}
}

void TestForViterbi(){
	InsertInitialValue("test.fasta");
	SetInitialForViterbi("abba");
	for(int i=0;i<state_size-1;i++){
		assert(hidenn_state_candidates_.at(i).at(0)==0);
	}
	assert(hidenn_state_percentages_.at(0).at(0)==log(0.2));
	assert(hidenn_state_percentages_.at(1).at(0)==log(0));
	assert(hidenn_state_percentages_.at(2).at(0)==log(0.06));
	CalculateAllStatesInViterbi("abba");
	//OutputHiddenStateForViterbi("abba");
}

void TestForForward(){
	InsertInitialValue("test.fasta");
	SetInitialForForward("abba");
	assert(forward_matrix_.at(0).at(0)>0.768&&forward_matrix_.at(0).at(0)<0.77);
	assert(forward_matrix_.at(1).at(0)==0);
	assert(forward_matrix_.at(2).at(0)>0.229&&forward_matrix_.at(2).at(0)<0.231);
	CalculateForForward("abba");
	assert(SumForward("abba")>0.99&&SumForward("abba")<1.01);
}

void LearnParameter(){
	sequence_={"GGGGATGTAGCTCAGATGGTAGAGCGCTCGCTTAGCATGCGAGAGGTACGGGGATCGATACCCCGCATCTCCA",
		"GATCGCATAGCGGAGTGGATATCGCGTTAGACTCCGAATCTAAAGGTCGTGGGTTCGAATCCCACTGCGATCG",
		"GGTCCCATGGTCTAGTGGTCAGGACATTGGACTCTGAATCCAGTAACCCGAGTTCAAATCTCGGTGGGACCT",
		"GTGGAGATGGCCGAGTTGGTCTAAGGCGCCAGATTAAGGTTCTTGTCCGAAAGGGCGTGGGTTCAAATCCCACTCTCCACA",
		"GGGGTGGTGGCGCAGTTGGCTAGCGCGTAGGTCTCATAGCTTCCGAGTAATCCTGAGGTCGAGAGTTCGAGCCTCTCTCACCCCA",
		"GGGCATTTGGTCTAGTGGTATGATTCTCGCTTTGGGTGCGAGAGGTCCCGAGTTCGATTCTCGGAATGACCC",
		"AGTCCCATAGCCTAGTTGGTCGAGCACAAGGTTTTCAATCTTGTGGTCGTGGGTTCGAGCCCTATTGGTGGTT",
		"CCGATCTTAACTCAATTGGTAGAGCAGAGGACCGTAGTGGGGTTGACCCAATAGATATCCTTAGGTCGTTAGTTTGAATCCAACAGGTCTAA",
		"GTCTAGGTGGTATAGTTGGTTATCACGCTAGTCTCACACACTAGAGGTCCCTAGTTCGAACCCAGGCTCAGATA",
		"GTCGATATGTCCGAGTGGTTAAGGAGACAGACTTGAAATCTGTTGGGCTTCGCCCGCGCAGGTTCGAACCCTGCTGTCGACG"
	};
	InsertInitialValue("initial.fasta");
	for(int j=0;j<800;j++){
		for(int i =0; i<sequence_.size();i++){
			string sequence=sequence_[i];
			SetInitialForForward(sequence);
			CalculateForForward(sequence);
			SetInitialForBackward(sequence);
			CalculateForBackward(sequence);
			Recalculate(sequence);
		}  
	}
}

void RunForKadai1(){
	string model="gagaguccuauacaaacuccaaaacacugagaccauacaaguaaaaccagucgagaaaauagucaacagcacccccguugcguauccugcggagaaugccucgcuaucuguccucacgauuggucuaaccgcucggcucaggcgugugggccugaaauccgggcagcaaacuacgguaaguuuucgcguaucaaaaauacauugaugaauguucgcuauuagccgggucgacguuuugauggugacuaggagcgaaagugauuuuuuugugagcggugucauagcaggggaucaucugaggugaacuauggacgggucagucgcccucauuggguuuguuacucagauacgugacacacguaaggucgcacggcaguagugauccacgagaaucggcacucuuacgagcaagucuauagcgacguggcuugcuauuaagacaguaauggacucggacagcuaguacgcgacuaaucauuuucgcacgccaucacacuuucaacccaggagcuuacaguguuccaggggccaggcuuggggcagccucuguggaagugcgagggcugcgcuaguguacauuagccgacccucagacgugaaauaaagagaugcugcugguugcacagugagcaacguucacucggcaacccgucggauacc";
	InsertInitialValue("initial.fasta");
	SetInitialForViterbi(model);
	CalculateAllStatesInViterbi(model);
	OutputHiddenStateForViterbi(model);
	cout<<endl;
}

void RunForKadai2(){
    
    sequence_={"GGGGATGTAGCTCAGATGGTAGAGCGCTCGCTTAGCATGCGAGAGGTACGGGGATCGATACCCCGCATCTCCA",
    "GATCGCATAGCGGAGTGGATATCGCGTTAGACTCCGAATCTAAAGGTCGTGGGTTCGAATCCCACTGCGATCG",
    "GGTCCCATGGTCTAGTGGTCAGGACATTGGACTCTGAATCCAGTAACCCGAGTTCAAATCTCGGTGGGACCT",
    "GTGGAGATGGCCGAGTTGGTCTAAGGCGCCAGATTAAGGTTCTTGTCCGAAAGGGCGTGGGTTCAAATCCCACTCTCCACA",
    "GGGGTGGTGGCGCAGTTGGCTAGCGCGTAGGTCTCATAGCTTCCGAGTAATCCTGAGGTCGAGAGTTCGAGCCTCTCTCACCCCA",
    "GGGCATTTGGTCTAGTGGTATGATTCTCGCTTTGGGTGCGAGAGGTCCCGAGTTCGATTCTCGGAATGACCC",
    "AGTCCCATAGCCTAGTTGGTCGAGCACAAGGTTTTCAATCTTGTGGTCGTGGGTTCGAGCCCTATTGGTGGTT",
    "CCGATCTTAACTCAATTGGTAGAGCAGAGGACCGTAGTGGGGTTGACCCAATAGATATCCTTAGGTCGTTAGTTTGAATCCAACAGGTCTAA",
    "GTCTAGGTGGTATAGTTGGTTATCACGCTAGTCTCACACACTAGAGGTCCCTAGTTCGAACCCAGGCTCAGATA",
    "GTCGATATGTCCGAGTGGTTAAGGAGACAGACTTGAAATCTGTTGGGCTTCGCCCGCGCAGGTTCGAACCCTGCTGTCGACG"
	};
    LearnParameter();
    for(int i =0; i<sequence_.size();i++){
    string sequence=sequence_[i];
    SetInitialForViterbi(sequence);
    CalculateAllStatesInViterbi(sequence);
    OutputHiddenStateForViterbi(sequence);
    cout<<endl;
    }
}
void RunForTestKadai2(){
    
    InsertInitialValue("test.fasta");
    string sequence="abba";
    SetInitialForForward(sequence);
	CalculateForForward(sequence);
	SetInitialForBackward(sequence);
	CalculateForBackward(sequence);
	Recalculate(sequence);
	assert(!state_transition_matrix[1][1]);
	assert(!state_transition_matrix[1][3]);
	assert(!state_transition_matrix[2][1]);
	assert(!state_transition_matrix[3][3]);
	assert(0.36597<state_transition_matrix[2][2]&&state_transition_matrix[2][2]<0.36598);
	assert(0.078349<state_transition_matrix[3][2]&&state_transition_matrix[3][2]<0.07835);
	assert(!output_probability_matrix[1][0]);//2→a
	assert(0.7391<output_probability_matrix[0][0]&&output_probability_matrix[0][0]<0.7392);
	assert(0.63<output_probability_matrix[2][1]&&output_probability_matrix[2][1]<0.631);//3→b
}
void RunForKadai3(){
    sequence_={"GGGGATGTAGCTCAGATGGTAGAGCGCTCGCTTAGCATGCGAGAGGTACGGGGATCGATACCCCGCATCTCCA",
    "GATCGCATAGCGGAGTGGATATCGCGTTAGACTCCGAATCTAAAGGTCGTGGGTTCGAATCCCACTGCGATCG",
    "GGTCCCATGGTCTAGTGGTCAGGACATTGGACTCTGAATCCAGTAACCCGAGTTCAAATCTCGGTGGGACCT",
    "GTGGAGATGGCCGAGTTGGTCTAAGGCGCCAGATTAAGGTTCTTGTCCGAAAGGGCGTGGGTTCAAATCCCACTCTCCACA",
    "GGGGTGGTGGCGCAGTTGGCTAGCGCGTAGGTCTCATAGCTTCCGAGTAATCCTGAGGTCGAGAGTTCGAGCCTCTCTCACCCCA",
    "GGGCATTTGGTCTAGTGGTATGATTCTCGCTTTGGGTGCGAGAGGTCCCGAGTTCGATTCTCGGAATGACCC",
    "AGTCCCATAGCCTAGTTGGTCGAGCACAAGGTTTTCAATCTTGTGGTCGTGGGTTCGAGCCCTATTGGTGGTT",
    "CCGATCTTAACTCAATTGGTAGAGCAGAGGACCGTAGTGGGGTTGACCCAATAGATATCCTTAGGTCGTTAGTTTGAATCCAACAGGTCTAA",
    "GTCTAGGTGGTATAGTTGGTTATCACGCTAGTCTCACACACTAGAGGTCCCTAGTTCGAACCCAGGCTCAGATA",
    "GTCGATATGTCCGAGTGGTTAAGGAGACAGACTTGAAATCTGTTGGGCTTCGCCCGCGCAGGTTCGAACCCTGCTGTCGACG"
    };
    LearnParameter();
    vector<vector<double>> initial_state_transition_matrix=state_transition_matrix;
    vector<vector<double>> initial_output_probability_matrix=output_probability_matrix;

    for(int i =0; i<sequence_.size()-1;i++){
	    state_transition_matrix=initial_state_transition_matrix;
	    output_probability_matrix=initial_output_probability_matrix;
	    string sequence1=sequence_[i];
	    SetInitialForForward(sequence1);
	    CalculateForForward(sequence1);
	    SetInitialForBackward(sequence1);
	    CalculateForBackward(sequence1);
	    //TODO Rewrite Recalcu funciyion for not using unnecessary function
	    Recalculate(sequence1);
	    vector<vector<double>> first_state_transition_=estimation_move_;
	    vector<vector<double>> first_output_probability_=estimation_output_;

	    for(int j=i+1;j<sequence_.size();j++){
		    state_transition_matrix=initial_state_transition_matrix;
		    output_probability_matrix=initial_output_probability_matrix;
		    string sequence2=sequence_[j];
		    SetInitialForForward(sequence2);
		    CalculateForForward(sequence2);
		    SetInitialForBackward(sequence2);
		    CalculateForBackward(sequence2);
		    Recalculate(sequence2);
		    double count_carnel=0;
		    for(int k=1;k<state_size;k++){
			    for(int l=1;l<state_size;l++){
				    count_carnel+=first_state_transition_[k][l]*estimation_move_[k][l];
			    }
		    }
		    for(int k=0;k<state_size;k++){
			    for(int l=0;l<alpha_size_;l++){
			    }
			    cout<<"sequence1="<<i<<"sequence2="<<j<<":"<<count_carnel<<endl;
		    }
	    }
}

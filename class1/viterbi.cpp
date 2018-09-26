
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <assert.h>
#include <math.h>
using namespace std;
int alpha_size_;
vector<char> alpha_;
int state_size;
vector<vector<int>> hidenn_state_candidates_;
vector<vector<double>> hidenn_state_percentages_;
vector<vector<double>> state_transition_matrix;
vector<vector<double>> output_probability_matrix;

double FindOutPutPercentage(char alpha,int state){
    for (int i=0;i<alpha_size_;i++) {
        if(alpha==alpha_[i]){
            return output_probability_matrix.at(state).at(i);
        }
    }
    return 0;
}


void SetInitialState(string model){
    hidenn_state_candidates_.clear();
    hidenn_state_percentages_.clear();
    vector<double> first_state_transition_matrix=*state_transition_matrix.begin();
    for (int i=0;i<state_size-1;i++) {
            vector<double> hidenn_state_percentages_line_;
            vector<int> hidenn_state_candidates_line_;
            hidenn_state_candidates_line_.push_back(0);
            double init_num=FindOutPutPercentage(model[0],i)*first_state_transition_matrix.at(i+1);
            hidenn_state_percentages_line_.push_back(init_num);
            hidenn_state_percentages_.push_back(hidenn_state_percentages_line_);
            hidenn_state_candidates_.push_back(hidenn_state_candidates_line_);
        }
}

void SetNotInitialState(string model){
    for (int model_number=1;model_number<model.size();model_number++){
        for (int after=1;after<state_size;after++) {
            int hidenn_state_max_candidate_=0;
            double hidenn_state_max_percentage_=0;
            for(int before=1;before<state_size;before++) {
                double state_percentage=
                hidenn_state_percentages_.at(before-1).at(model_number-1)*state_transition_matrix.at(before).at(after)*FindOutPutPercentage(model[model_number],after-1);
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

void OutputHiddenState(string model){
    vector<int> hidden_state_;
    double hidden_state_max_percentage_=0;
    int hidden_state_max_candidate_=0;
    for(int i=1;i<state_size;i++){
        double hidden_state_max_percentage_candidate_=hidenn_state_percentages_.at(i-1).back();
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

void InsertInitialValueInLog(string file){
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
void TestForStates(){
    InsertInitialValueInLog("test.fasta");
    SetInitialState("abba");
    for(int i=0;i<state_size-1;i++){
        assert(hidenn_state_candidates_.at(i).at(0)==0);
    }
    assert(hidenn_state_percentages_.at(0).at(0)==0.2);
    assert(hidenn_state_percentages_.at(1).at(0)==0);
    assert(hidenn_state_percentages_.at(2).at(0)==0.06);
    SetNotInitialState("abba");
    OutputHiddenState("abba");
    
}

int main(void){
    //TestForStates();
    string model="gagaguccuauacaaacuccaaaacacugagaccauacaaguaaaaccagucgagaaaauagucaacagcacccccguugcguauccugcggagaaugccucgcuaucuguccucacgauuggucuaaccgcucggcucaggcgugugggccugaaauccgggcagcaaacuacgguaaguuuucgcguaucaaaaauacauugaugaauguucgcuauuagccgggucgacguuuugauggugacuaggagcgaaagugauuuuuuugugagcggugucauagcaggggaucaucugaggugaacuauggacgggucagucgcccucauuggguuuguuacucagauacgugacacacguaaggucgcacggcaguagugauccacgagaaucggcacucuuacgagcaagucuauagcgacguggcuugcuauuaagacaguaauggacucggacagcua";
    //guacgcgacuaaucauuuucgcacgccaucacacuuucaacccaggagcuuacaguguuccaggggccaggcuuggggcagccucuguggaagugcgagggcugcgcuaguguacauuagccgacccucagacgugaaauaaagagaugcugcugguugcacagugagcaacguucacucggcaacccgucggauacc ";
    InsertInitialValueInLog("initial.fasta");
    SetInitialState(model);
    SetNotInitialState(model);
    OutputHiddenState(model);
}

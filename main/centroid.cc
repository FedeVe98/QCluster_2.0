#include <string.h>
#include <math.h>
#include<iostream>
#include <unordered_map>
using namespace std;

// mean centroid
void mean_centroid(int N, int row_length,  unordered_map<string, double>** freq, const int num_nt,
                   double** freq_1, unordered_map<string, double>* centroid, unordered_map<string, double>* centroid_tilde,
                   unordered_map<string, double>**quality,  unordered_map<string, double>*expected_qual, double **quality_1,
                   unordered_map<string, double>*expected_freq)
{
	for(int n=0; n<N; n++){		
		for(auto iter = quality[n]->begin(); iter != quality[n]->end(); ++iter){
			if(centroid->find(iter->first) == centroid->end())
				centroid->insert(make_pair(iter->first, iter->second));
			else
				centroid->at(iter->first) += iter->second;
		}
	}
	for(auto iter = centroid->begin(); iter != centroid->end(); ++iter)
		centroid->at(iter->first) /= N;
		
	return;
}


// d2 centroid
void d2_centroid(int N, int row_length,  unordered_map<string, double>** freq, const int num_nt,
                 double** freq_1, unordered_map<string, double>* centroid, unordered_map<string, double>* centroid_tilde,
                 unordered_map<string, double>**quality,  unordered_map<string, double>*expected_qual, double **quality_1,
                 unordered_map<string, double>*expected_freq)
{
	for(int n=0; n<N; ++n){
		for(auto iter = quality[n]->begin(); iter != quality[n]->end(); ++iter){
			if(centroid->find(iter->first) == centroid->end())
				centroid->insert(make_pair(iter->first, iter->second));
			else
				centroid->at(iter->first) += iter->second;
		}
	}
	// now normalize the entries
	double total_sq_count = 0;
	for(auto iter = centroid->begin(); iter != centroid->end(); ++iter){
		total_sq_count += iter->second * iter->second;
	}
	double norm = sqrt(total_sq_count);
	for(auto iter = centroid->begin(); iter != centroid->end(); ++iter)
		centroid->at(iter->first) /= norm;
	return;
}


// KL centroid --- operates on raw counts which do not have to be normalized
void kl_centroid(int N, int row_length,  unordered_map<string, double>** freq, const int num_nt,
                 double** freq_1, unordered_map<string, double>* centroid, unordered_map<string, double>* centroid_tilde,
                 unordered_map<string, double>**quality,  unordered_map<string, double>*expected_qual, double **quality_1,
                 unordered_map<string, double>*expected_freq)
{
	for(int n=0; n<N; ++n){
		for(auto iter = quality[n]->begin(); iter != quality[n]->end(); ++iter){
			if(centroid->find(iter->first) == centroid->end())
				centroid->insert(make_pair(iter->first, iter->second));
			else
				centroid->at(iter->first) += iter->second;
		}
	}
	// now normalize the entries
	double total_count = 0;
	for(auto iter = centroid->begin(); iter != centroid->end(); ++iter)
		total_count += iter->second;
	for(auto iter = centroid->begin(); iter != centroid->end(); ++iter)
		centroid->at(iter->first) /= total_count;
	return;
}


// MM centroid --- evaluates frequencies of single nucleotides
// and coumputes frequencies of words using zero order Markov model;
// i.e., as product of single nucleotide frequencies:
// f_{w1...w_k} = f_w1 f_w2 ... f_wk
void mm_centroid(int N, int row_length,  unordered_map<string, double>** freq, const int num_nt,
                 double** freq_1, unordered_map<string, double>* centroid, unordered_map<string, double>* centroid_tilde,
                 unordered_map<string, double>**quality,  unordered_map<string, double>*expected_qual, double **quality_1,
                 unordered_map<string, double>*expected_freq)
{
	//If p_method is global (P1G or P2G) we don't need to calculate centroids
	if (expected_freq != NULL) return;
		
	//P2L: It uses a markovian model to compute the expected frequency of words
	double nt_freq[num_nt];
	memset(nt_freq, 0, num_nt * sizeof(nt_freq[0]));
	double S=0;
	// compute frequencies of individual nucleotides
	for (int n=0; n<N; n++){
		for (int l=0; l<num_nt; l++){
			//Number of occurrencies of LETTERS in the cluster:
			nt_freq[l] += freq_1[n][l]; 
		}
	}
	for (int l=0; l<num_nt; l++)
		S += nt_freq[l];	//S = number of letter in the cluster
	for (int l=0; l<num_nt; l++)
		nt_freq[l] /= S;	//Average frequency of letters in the cluster

	// compute frequencies of words using Markov model
	for(int l=0; l<row_length; l++){
		double p = 1;
		for(int M=1; M<row_length; M*=num_nt){
			int digit = l/M % num_nt;
			p *= nt_freq[digit];
		}
		//centroid[l] = p;
	}
	return;
}


// d2* centroid
void d2ast_centroid(int N, int row_length,  unordered_map<string, double>** freq, const int num_nt,
                    double** freq_1, unordered_map<string, double>* centroid, unordered_map<string, double>* centroid_tilde,
                    unordered_map<string, double>**quality,  unordered_map<string, double>*expected_qual, double **quality_1,
                    unordered_map<string, double>*expected_freq)
{
	// compute frequencies from MM
	mm_centroid(N, row_length, freq, num_nt, freq_1, centroid, NULL, 
				quality, expected_qual, quality_1, expected_freq);
	double total_count = 0;
	for(int n=0; n<N; ++n){
		for(auto iter = quality[n]->begin(); iter != quality[n]->end(); ++iter){
			if(centroid->find(iter->first) == centroid->end())
				centroid_tilde->insert(make_pair(iter->first, iter->second));
			else
			{
				centroid_tilde->at(iter->first) += quality[n]->operator[](iter->first);
				total_count += freq[n]->operator[](iter->first);
			}
		}
		
	}
	
	//Determines what data to use: global(expected_freq) or local(centroid)
	unordered_map<string, double> * exp_freq;
	if (expected_freq==NULL) exp_freq = centroid;
	else exp_freq = expected_freq;
	
	double S = 0;
	for(auto iter = centroid_tilde->begin(); iter != centroid->end(); ++iter){
		centroid_tilde->at(iter->first) -= total_count * exp_freq->operator[](iter->first) * expected_qual->operator[](iter->first);
		centroid_tilde->at(iter->first) /= sqrt(exp_freq->operator[](iter->first) * expected_qual->operator[](iter->first));
		S += iter->second * iter->second;
	}
	S = sqrt(S);
	for(auto iter = centroid_tilde->begin(); iter != centroid->end(); ++iter)
		centroid_tilde->at(iter->first) /= S;
	return;
}

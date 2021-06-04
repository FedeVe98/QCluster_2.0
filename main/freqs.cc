#include <math.h>
#include "freqs.hh"
#include <iostream>
#include <cstdlib>
#include <string.h>
#include <unordered_map>
using namespace std;

//UTILITY FUNCTION
static bool is_valid_nt(char c)
{
	c = toupper(c);
	if (c == 'A') {
		return true;
	} 
	if (c == 'C') {
		return true;
	} 
	if (c == 'G') {
		return true;
	} 
	if (c == 'T') {
		return true;
	}
	return false;
}


static int nt2int(char c)
{
	switch (c) {
		case 'A':
		case 'a': 
			return 0;
		case 'C':
		case 'c': 
			return 1;
		case 'G':
		case 'g': 
			return 2;
		case 'T':
		case 't': 
			return 3;
	}
	// should never get here
	return -1;
}

static char int2nt(int c)
{
	switch (c) {
		case 0: 
			return 'A';
		case 1: 
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
	// should never get here
	return -1;
}

//NORMALIZE FUNCTIONS
//Divide each element of the vector by the norm1 of the vector itself
void normalize_row(unordered_map<string, double>* fv)
{
	double total = 0;
	for (auto iter = fv->begin(); iter != fv->end(); ++iter){
		total += iter->second;	// ||fv||1
	}
	for (auto iter = fv->begin(); iter != fv->end(); ++iter){
		fv->operator[](iter->first) /= total;
	}
	return;
}

//Divide each element of "qual" vector by the norm1 of "freq" vector
//Freq vector MUST NOT be normalized
void normalize_row_against(unordered_map<string, double>* freq, unordered_map<string, double>* qual)
{
	double total = 0;
	for (auto iter = freq->begin(); iter != freq->end(); ++iter){
		total += iter->second;
	}
	for (auto iter = qual->begin(); iter != qual->end(); ++iter){
		qual->operator[](iter->first) /= total;
	}
	return;
}

//Project Freq vector on the hyperplane x1 + x2 + x3... -1 = 0
//Freq and qual vector MUST be normalized (with normalize_row_against)
void normalize_row_projection(unordered_map<string, double>* freq, unordered_map<string, double>* qual)
{
	double total = 0, total_qual = 0;
	for (auto iter = freq->begin(); iter != freq->end(); ++iter){
		total += iter->second;
		total_qual += qual->operator[](iter->first);
	}
	double partial = (total - total_qual) / freq->size();
	for (auto iter = qual->begin(); iter != qual->end(); ++iter){
		qual->operator[](iter->first) += partial;
	}
	return;
}

//FREQUENCY FUNCTIONS
//Note: freq_1 MUST NOT be normalized
//nucleotides of kmers
void fill_overlap_count_vector_1(string seq, string seq_qual, int K,
                                double* freq_1, double* quality_vector, int normalize, double pc,
                                double *avg_quality_1, bool redistribute)
{
	int L = seq.length();
	int *kmer = new int[K];
	int *qual = new int[K];
	double *prob = new double[K];
	double *freq_bases = new double[NUM_NT];
	double readqual=1;
	for(int i=0; i<L; ++i){
		seq[i] = toupper(seq[i]);
	}
	//If K==1, freq_1 is NULL and freq_bases is initialized to 1/4 for each base
	if (K==1){
		for (int i=0; i<NUM_NT; i++) freq_bases[i] = 1.0/4.0;
	}
	bool valid_kmer = true;
	int index_kmer = 0;
	string kmer_string;
	for (int i=0; i<L-K+1; i++){
		valid_kmer = true;
		index_kmer = 0;
		readqual = 1;
		//Fill the kmer, qual and prob vectors
		for (int j=0; j<K; j++){
			if (!is_valid_nt(seq[i+j])){
					valid_kmer = false;
					break;
				}
			kmer[j] = nt2int(seq[i+j]);
			kmer_string += seq[i+j];
			qual[j] = int(seq_qual[i+j])-base;
			prob[j] = (1.0 - pow(10.0,-(qual[j])/10.0));
			readqual *= prob[j];
            index_kmer = NUM_NT*index_kmer + kmer[j];
		}
		if (!valid_kmer) continue;
		freq_1[index_kmer] += 1;
		quality_vector[index_kmer] += readqual;
		if (avg_quality_1 != NULL) 
			avg_quality_1[kmer[0]] += qual[0];
	}

	/*** PRINTING FREQ_1 TEST ***/
	/*int count = 0;
	for(int i = 0; i < L; i++)
	{
		cout << freq_1[i] << endl;
		count++;

		if(count>10)	continue;
	}*/

	delete[] kmer;
	delete[] qual;
	delete[] prob;
	if (K==1) delete[] freq_bases;
	return;
}

//kmers of sequences
void fill_overlap_count_vector(string seq, string seq_qual, int K,
						unordered_map<string, double>* freq, unordered_map<string, double>* quality_vector,
						int normalize, double pc, unordered_map<string, double*> *avg_quality,
						double* freq_1, bool redistribute)
{
	int L = seq.length();
	int *kmer = new int[K];
	int *qual = new int[K];
	double *prob = new double[K];
	double *freq_bases = new double[NUM_NT];
	double readqual=1;
	for(int i=0; i<L; ++i){
		seq[i] = toupper(seq[i]);
	}
	bool valid_kmer = true;
	int index_kmer = 0;
	string kmer_string;

	//If K==1, freq_1 is NULL and freq_bases is initialized to 1/4 for each base
	if (K==1){
		for (int i=0; i<NUM_NT; i++) freq_bases[i] = 1.0/4.0;
	}//Otherwise, freq_bases is the normalization of freq_1
	else{
		double somma =0;
		for (int i=0; i<NUM_NT; i++) somma += freq_1[i];
		for (int i=0; i<NUM_NT; i++) freq_bases[i] = freq_1[i]/somma;
	}
	for (int i=0; i<L-K+1; i++){
		valid_kmer = true;
		readqual = 1;
		kmer_string = "";
		//Fill the kmer, qual and prob vectors
		for (int j=0; j<K; j++){
			if (!is_valid_nt(seq[i+j])){
					valid_kmer = false;
					break;
				}
			kmer[j] = nt2int(seq[i+j]);
			kmer_string += seq[i+j];
			qual[j] = int(seq_qual[i+j])-base;
			prob[j] = (1.0 - pow(10.0,-(qual[j])/10.0));
			readqual *= prob[j];
		}
		if (!valid_kmer) continue;
		//insert into freq data structure
		if(freq->find(kmer_string) == freq->end())
			freq->insert(make_pair(kmer_string, pc));
		else
			freq->operator[](kmer_string) += 1;
		quality_vector->operator[](kmer_string) += readqual;

		if (avg_quality != NULL)
		{
			if(avg_quality->find(kmer_string) == avg_quality->end())
			{
				double temp_array[NUM_NT] = {};
				avg_quality->insert(make_pair(kmer_string, &temp_array[0]));

			}

			for (int j=0; j<K; j++) avg_quality->operator[](kmer_string)[kmer[j]] += qual[j];
		}

		if (!redistribute) continue;

		//Redistributing quality: for each letter in the kmer...
		kmer_string="";
		for (int j=0; j<K; j++){
			int original_letter = kmer[j];
			double original_prob = prob[j];
			//...we try to replace it with another
			for (int letter=0; letter<NUM_NT; letter++){
				//we must redistribute also to the same letter
				kmer[j] = letter;
				//Equally subdivided probability
				//Probability subidivided basing on the letter frequency
				prob[j] = (1.0-original_prob) * double(freq_bases[letter]);
				//Recalculation of the new kmer probability
				readqual = 1;
				index_kmer = 0;
				for (int f=0; f<K; f++){
					readqual *= prob[j];
					char tmp = int2nt(kmer[j]);
                    kmer_string += tmp;
				}
				//Add a little probability to that kmer
				quality_vector->operator[](kmer_string) += readqual;
			}
			kmer[j] = original_letter;
			prob[j] = original_prob;
		}
	}


	int N = freq->size();
	switch(normalize){
		case 1:	
			normalize_row(freq);
			normalize_row(quality_vector);
			break;
		case 2:
			normalize_row_against(freq, quality_vector);
			normalize_row(freq);
			break;
		case 3:
			normalize_row_against(freq, quality_vector);
			normalize_row(freq);
			normalize_row_projection(freq, quality_vector);
			break;
		default:
			break;
	}

	/*** PRINTING FREQ_1 TEST ***/
	/*int count = 0;
	for(auto iter = freq->begin(); iter != freq->end(); ++iter)
	{
		cout << iter->first << "  " << iter->second << endl;
		count++;

		if(count>10)	continue;
	}*/
	
	delete[] kmer;
	delete[] qual;
	delete[] prob;
	if (K==1) delete[] freq_bases;
	return;
}




/**Function responsible for the calculation of the quality expected value for 
 * each k-word: the expected value is computed basing on the full dataset.
 * Refer to the documentation for the description of the possible methods. 
 * Default method: E1 */
void calculate_quality_expected_value(int method, int N, int K, int L,
                                      unordered_map<string, double> **freq, unordered_map<string, double> **quality,
                                      double **freq_1, double* avg_quality_1, unordered_map<string,
                                      double*> *avg_quality, double *expected_qual_1, unordered_map<string, double> *expected_qual){
	switch (method){
	case 2:
		{
			int *times_to_use = new int[NUM_NT];
			for(auto iter = avg_quality->begin(); iter != avg_quality->end(); ++iter)
			{
				for(int j=0; j<NUM_NT; j++) times_to_use[j] = 0;

				int index;
				for(int k=0; k<K; k++)
				{
					index = nt2int(iter->first[k]);
					times_to_use[index] += 1;
				}

				for(int j=0; j < NUM_NT; j++)
				{
					if(times_to_use[j] != 0)
					{
						avg_quality->operator[](iter->first)[j] /= times_to_use[j];
					}
				}
			}

			delete[] times_to_use;

			for(auto iter = avg_quality->begin(); iter != avg_quality->end(); ++iter)
			{
				expected_qual->insert(make_pair(iter->first, 1));

				double num_occ = 0;

				for(int i=0; i<N; i++)
				{
					if(freq[i]->find(iter->first) != freq[i]->end())
						num_occ += freq[i]->operator[](iter->first);
				}

				for(int j=0; j < NUM_NT; j++)
				{
					if(avg_quality->operator[](iter->first)[j] != 0)
					{
						expected_qual->operator[](iter->first) *= 1 -
						pow(10.0, -(avg_quality->operator[](iter->first)[j])
										/(10.0*num_occ));
					}
				}
			}
		break;
		}
	case 3:
		{
		for(int i=0; i<NUM_NT; i++) {
			int tot = 0;
			for (int j=0; j<N; j++) tot += freq_1[j][i];
			avg_quality_1[i] = 1-pow(10.0,-avg_quality_1[i]/(10.0*tot));
		}
		int divisore = 1;
		
		for (int i=0; i<L; i++) expected_qual_1[i] = 1;
		for(int k=0; k<K; k++){
			for (int i=0; i<L; i++){
				expected_qual_1[i] *= avg_quality_1[(i/divisore)%(NUM_NT)];
			}
			divisore *= NUM_NT;
		}
		break;
		}
	case 1: //fall
		
	default:
		{
		unordered_map<string, double> *denominator = new unordered_map<string, double>();
		for(int i = 0; i < N; i++)
		{
			for(auto iter = freq[i]->begin(); iter != freq[i]->end(); ++iter)
			{
				if(expected_qual->find(iter->first) == expected_qual->end())
					expected_qual->insert(make_pair(iter->first, 0));

				expected_qual->operator[](iter->first) += quality[i]->operator[](iter->first);
				
				if(denominator->find(iter->first) == denominator->end())
					denominator->insert(make_pair(iter->first, 0));
				
				denominator->operator[](iter->first) += iter->second;
			}
		}

		for(auto iter = denominator->begin(); iter != denominator->end(); ++iter)
			expected_qual->operator[](iter->first) /= iter->second;


		break;	
		}
	}
}

//expected frequency of nucleotides
void expected_frequency_p2global(int N, int K, int L, double *expected_freq, double **freq_1)
{
	double nt_freq[NUM_NT];
	memset(nt_freq, 0, NUM_NT * sizeof(nt_freq[0]));
	double S=0;
	// compute frequencies of individual nucleotides
	for (int n=0; n<N; n++){
		for (int l=0; l<NUM_NT; l++){
			nt_freq[l] += freq_1[n][l]; //= count of letter l on the dataset
		}
	}
	for (int l=0; l<NUM_NT; l++)
		S += nt_freq[l];	//S = total number of letters in the dataset

	for (int l=0; l<NUM_NT; l++){
		nt_freq[l] /= S;	//Average frequency of each letter in the dataset
	}
	// compute frequencies of words using Markov model
	for(int l=0; l<L; l++){
		double p = 1;
		for(int M=1; M<L; M*=NUM_NT){
			int digit = l/M % NUM_NT;
			p *= nt_freq[digit];
		}
		expected_freq[l] = p;
	}
	return;
}

//expected frequency of kmers
void expected_frequency_p1global(int N, int K, int L, unordered_map<string, double> *expected_freq,
								 unordered_map<string, double> **freq)
{
	int tot_words = 0;
	for(int i = 0; i < N; i++)
	{
		for(auto iter = freq[i]->begin(); iter != freq[i]->end(); ++iter)
		{
			if(expected_freq->find(iter->first) == expected_freq->end())
				expected_freq->insert(make_pair(iter->first, 0));

			expected_freq->operator[](iter->first) += iter->second;
			
			tot_words += expected_freq->operator[](iter->first);
		}
	}
	for(auto iter = expected_freq->begin(); iter != expected_freq->end(); ++iter)
		expected_freq->operator[](iter->first) /= tot_words;
	return;
}

// normalize frequency matrix to make its columns univariant
void normalize_freq_matrix(unordered_map<string, double> **freq, unordered_map<string, double> **qual, int N, int row_length)
{
	double tmp, sum_freq, sum_freq_square, V_freq;
    double sum_qual, sum_qual_square, V_qual;
    for(int i = 0; i < N; i++)
            {
                sum_freq_square=0;
				sum_freq=0;
				sum_qual_square = 0;
				sum_qual = 0;
				for(auto iter = freq[i]->begin(); iter != freq[i]->end(); ++iter)
				{
					tmp = iter->second;
					sum_freq_square += tmp*tmp;
					sum_freq += tmp;
					tmp = qual[i]->operator[](iter->first);
					sum_qual_square += tmp*tmp;
					sum_qual += tmp;
				}
				//Variance
				V_freq = sqrt( sum_freq_square/N - (sum_freq*sum_freq)/(N*N) );
				V_qual = sqrt( sum_qual_square/N - (sum_qual*sum_qual)/(N*N) );
				if (V_freq>0 && V_qual>0){
					for(auto iter = freq[i]->begin(); iter != freq[i]->end(); ++iter){
						iter->second /= V_freq;
						iter->second /= V_qual;
					}
				}
            }
    return;
}

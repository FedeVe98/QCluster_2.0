#include <string.h>
#include <math.h>
#include<iostream>
#include <unordered_map>
using namespace std;



double euclidean_distance(unordered_map<string, double> *centroid, unordered_map<string, double>*y,
                          int L, unordered_map<string, double>* xt, unordered_map<string, double>* quality,
                          unordered_map<string, double>*expected_qual, unordered_map<string, double> *expected_freq)
{
	double S = 0;
	double z;
    for(auto iter = quality->begin(); iter != quality->end(); ++iter){
            z = centroid->operator[](iter->first) - iter->second;
            S += z*z;
    }
	return S;
}


double kl_distance(unordered_map<string, double> *centroid, unordered_map<string, double>*p, int L,
                   unordered_map<string, double>* centroid_tilde, unordered_map<string, double>* quality,
                   unordered_map<string, double>*expected_qual, unordered_map<string, double> *expected_freq)
{
	double total_count = 0;
    for(auto iter = quality->begin(); iter != quality->end(); ++iter){
            total_count += iter->second;
    }
	double S = 0;
    for(auto iter = quality->begin(); iter != quality->end(); ++iter){
            S += iter->second * log(iter->second/(total_count * centroid->operator[](iter->first)));
    }
	return S;
}


double symkl_distance(unordered_map<string, double> *centroid, unordered_map<string, double>*p, int L,
                      unordered_map<string, double>* centroid_tilde, unordered_map<string, double>* quality,
                      unordered_map<string, double>*expected_qual, unordered_map<string, double> *expected_freq)
{
	double total_count = 0;
    for(auto iter = quality->begin(); iter != quality->end(); ++iter) {
        total_count += iter->second;
    }
	double S = 0;
    for(auto iter = quality->begin(); iter != quality->end(); ++iter){
            S += (iter->second - total_count*centroid->operator[](iter->first) * log(iter->second/(total_count * centroid->operator[](iter->first))));
    }
	return S;
}


double d2_distance(unordered_map<string, double> *centroid, unordered_map<string, double>*p, int L,
                   unordered_map<string, double>* centroid_tilde, unordered_map<string, double>* quality,
                   unordered_map<string, double>*expected_qual, unordered_map<string, double> *expected_freq)
{
	double norm = 0;
	double S = 0;
    for(auto iter = quality->begin(); iter != quality->end(); ++iter){
        norm += iter->second * iter->second;
    }
	norm = sqrt(norm);
    for(auto iter = quality->begin(); iter != quality->end(); ++iter){
            S += (centroid->operator[](iter->first) * iter->second )/ norm;
    }
	return 1 - S;
}


double chi2_distance(unordered_map<string, double> *centroid, unordered_map<string, double>*p, int L,
                     unordered_map<string, double>* centroid_tilde, unordered_map<string, double>* quality,
                     unordered_map<string, double>*expected_qual, unordered_map<string, double> *expected_freq)
{
	double chi2 = 0;
	double total_count = 0;
    for(auto iter = quality->begin(); iter != quality->end(); ++iter){
            total_count += iter->second;
    }
    for(auto iter = quality->begin(); iter != quality->end(); ++iter){
            double exp_count = centroid->operator[](iter->first) * total_count;
            chi2 += (iter->second - exp_count) * (iter->second - exp_count) / exp_count;
    }
	//return (chi2 - (L-3) * log(chi2))/2;
	return chi2;
}

double d2ast_distance(unordered_map<string, double> *centroid, unordered_map<string, double>*p, int L,
                      unordered_map<string, double>* centroid_tilde, unordered_map<string, double>* quality,
                      unordered_map<string, double>*expected_qual, unordered_map<string, double> *expected_freq)
{
    //Determines what data to use: global(expected_freq) or local(centroid)
    unordered_map<string, double> * exp_freq;
    double total_count = 0;
    unordered_map<string, double>* x = new unordered_map<string, double>();
    //copy of quality map
    for (auto iter = quality->begin(); iter != quality->end(); ++iter)
        x->insert(make_pair(iter->first, iter->second));

    double S = 0;
    if (expected_freq==NULL) 
        exp_freq = centroid;
    else 
        exp_freq = expected_freq;
    for(auto iter = p->begin(); iter != p->end(); ++iter){
        total_count += iter->second;
    }
    for(int i=0; i<L; i++)
    {
        for(auto iter = expected_qual[i].begin(); iter != expected_qual[i].end(); ++iter){
            x->at(iter->first) -= total_count * exp_freq->operator[](iter->first) * iter->second;
            x->at(iter->first) /= sqrt(exp_freq->operator[](iter->first) * iter->second );
            S += x->operator[](iter->first) * x->operator[](iter->first);
        }
    }
    
    S = sqrt(S);
    for(auto iter = x->begin(); iter != x->end(); ++iter){
        x->at(iter->first) /= S;
    }
    S = 0;
    for(auto iter = x->begin(); iter != x->end(); ++iter){
        S += x->operator[](iter->first) * centroid_tilde->operator[](iter->first);
    }
    return (1-S)/2;
}


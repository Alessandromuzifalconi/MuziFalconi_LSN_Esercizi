#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "statistics.h"

using namespace std;

double error (double* ave, double* ave2, int i) {
	return sqrt((ave2[i]-pow(ave[i],2))/i);
}

double mean (int N, double* x) {
	double sum;
	sum = 0;
	for(int i=0; i<N; i++){
    	sum = sum + x[i]; 
   }
   return sum/N;
}

double Xi2 (int N, double* x, double exp) {
	double sum = 0;
	for(int i = 0; i<N; i++){
    	sum = sum + pow((x[i]-exp),2)/exp; 
   }
   return sum;
}

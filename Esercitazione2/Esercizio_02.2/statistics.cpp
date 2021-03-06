#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "statistics.h"

using namespace std;

double error (double ave, double ave2, int i) {
	if (i == 0||ave == ave2)
		return 0;
	else
		return sqrt((ave2-pow(ave,2))/double(i));
}

double mean (int N, double* x) {
	double sum;
	sum = 0;
	for(int i=0; i<N; i++){
    	sum = sum + x[i]; 
   }
   return sum/N;
}

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "statistics.h"

using namespace std;

double mean (int N, double x[]) {
	double sum;
	sum = 0;
	for(int i=0; i<N; i++){
    	sum = sum + x[i]; 
   }
   return sum/N;
}

double StdDev(int N, double x[]) {
    double m = mean(N,x);
    double sum2 = 0;
    for (int k = 0; k < N; k++)
        sum2 += pow(x[k]-m,2);
    sum2 = sum2/(N-1);
    return sqrt(sum2);
}

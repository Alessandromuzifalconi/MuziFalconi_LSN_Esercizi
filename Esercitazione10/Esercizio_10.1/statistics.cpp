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


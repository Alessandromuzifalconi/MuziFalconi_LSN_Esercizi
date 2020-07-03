#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistics.h"
#include "functions.h" 

using namespace std;

int main (int argc, char *argv[]){
	
	ofstream uniform_results, importance_results;
	
	ifstream Primes("Primes"); 
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
		input >> property;
		if( property == "RANDOMSEED" ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			rnd.SetRandom(seed,p1,p2);
		}
		}
	input.close();
   	} 
   	else cerr << "PROBLEM: Unable to open seed.in" << endl;
   	
	
	for (int i = 0; i < N; i++) { 					//Cycle over different blocks
		sum = 0;						
		sum_imp = 0;
		for (int j = 0; j < L; j++) {				//Cycle inside a block
			x = rnd.Rannyu();						//Generate a uniform random number
			sum += cos(M_PI/2*x)*M_PI/2;			//Accumulate the generated random numbers
			x = importance_sampl();					//Generate a random number distributed as p(x) = 2(1-x)
			sum_imp += M_PI/2*cos(M_PI/2*x)/(2*(1-x));							
		}
		ave[i] = sum/L;
		ave2[i] = pow(ave[i],2);
		ave_imp[i] = sum_imp/L;
		ave2_imp[i] = pow(ave_imp[i],2);
	}
	
	uniform_results.open("uniform_results.dat");
   	importance_results.open("imp_results.dat");
	
	for (int i = 0; i < N; i++) {    				//Blocking method
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
			sum_prog_imp[i] += ave_imp[j];
			sum2_prog_imp[i] += ave2_imp[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		sum_prog_imp[i] = sum_prog_imp[i]/(i+1);
		sum2_prog_imp[i] = sum2_prog_imp[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		err_prog_imp[i] = error(sum_prog_imp, sum2_prog_imp, i);
		
		uniform_results << sum_prog [i] << "\t" << err_prog[i] << endl; 			
		importance_results << sum_prog_imp [i] << "\t" << err_prog_imp[i] << endl; 	
	}
	
	uniform_results.close();
	importance_results.close();
	rnd.SaveSeed();
return 0;
}

double importance_sampl () { 				//Returns a random number distributed as p(x) = 2(1-x)
	double x = rnd.Rannyu();
	return 1-sqrt(1-x);
}

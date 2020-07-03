#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistics.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

	ofstream Mean_results, Stdev_results, Xi2_results;
	
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
	
		
	for (int i = 0; i < N; i++) { 		//Cycle over different blocks
		sum = 0;						
		sum_var = 0;
		for (int j = 0; j < L; j++) {	//Cycle inside a block
			x = rnd.Rannyu();			//Generate a random number
			sum += x;					//Accumulate the generated random numbers
			sum_var += pow(x-0.5,2);	//Accumulate the squared distances from the expectation value
		}
		ave[i] = sum/L;
		ave2[i] = pow(ave[i],2);
		ave_var[i] = sum_var/L;
		ave2_var[i] = pow(ave_var[i],2);
	}
	
	Mean_results.open("mean_results.dat");
   	Stdev_results.open("stdev_results.dat");
	
	for (int i = 0; i < N; i++) {    		//Blocking method for computing the mean value and the variance with respective uncertainties
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
			sum_prog_var[i] += ave_var[j];
			sum2_prog_var[i] += ave2_var[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		sum_prog_var[i] = sum_prog_var[i]/(i+1);
		sum2_prog_var[i] = sum2_prog_var[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		err_prog_var[i] = error(sum_prog_var, sum2_prog_var, i);
		
		Mean_results << sum_prog [i] << "\t" << err_prog[i] << endl; 			//Write results on a file
		Stdev_results << sum_prog_var [i] << "\t" << err_prog_var[i] << endl; 	
	}
	
	Mean_results.close();
	Stdev_results.close();
	cout << "Estimation of the epectation value and variance done, Xi2 starting..." << endl;
	//Compute Xi2	   	
   	for (int i = 0; i<m; i++)	//Initialize the counters to 0
   		count[i] = 0;
   		
   	for (int k = 0; k < N; k++) { 			//Cycle over blocks
   		for (int j = 0; j < m; j++) { 		//Cycle over the sub-intervals j
   			counter = 0; 					//Initialize the counter for the considered sub-interval to 0
   			for (int i = 0; i < n; i++) { 	//Generate n pseudo random numbers
   				x = rnd.Rannyu();
   				if (x > 1.0/m*j && x < 1.0/m*(j+1)) { 	//Check if the generated number is inside the right sub-interval and increase the counter if it is
   					counter = counter + 1;
				}
   			}
   			count[j] = counter; 						//Save the value of the counter for the considered sub-interval in the counter vector
   		}
   		Xi2_block[k] = Xi2(m, count, (double)n/m); 		//Xi2 for the considered block
   		Xi2_block2[k] = pow(Xi2_block[k], 2);			//Squared Xi2 for the considered block
	}
	
	for (int i = 0; i < N; i++) {						//Initialize all the sum variables to 0
		sum_prog[i] = 0;
		sum2_prog[i] = 0;
		err_prog[i] = 0;
	}
	
	Xi2_results.open("Xi_results.dat");
	 	  	
	 for (int i = 0; i < N; i++) {						//Blocking method for computing the Xi2 and its uncertainty 
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += Xi2_block[j];
			sum2_prog[i] +=  Xi2_block2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);	
		err_prog[i] = error(sum_prog, sum2_prog, i);
		Xi2_results << sum_prog [i] << "\t" << err_prog[i] << endl; 	//Write results on a file
	}	  	
	 	  
   	Xi2_results.close();
	rnd.SaveSeed();
return 0;
}

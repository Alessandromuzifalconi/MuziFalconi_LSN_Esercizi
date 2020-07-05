#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistics.h"
#include "functions.h"

using namespace std;


int main (int argc, char *argv[]){

	ofstream output_final, output_path;
	
	Random_Initializer();
	
	//Compute the option prices directly
	for (int i = 0; i < steps; i++) {
		for (int j = 0; j < L; j++) {
			ST = S0*exp((r-pow(sigma,2)/2)*T + sigma*rnd.Gauss(0.,sqrt(T))); 	//Directly evaluate the price S(T)
        	//Choose whether to evaluate the call or put option: 
        	if (0 > (ST-K)) {													//If the strike price is more than the actual price the holder sells the asset (put-option)
        	    C[i] = 0;
        		P[i] = exp(-r*T) *(K-ST);
       		} else {															//Otherwise the holder buys
            	C[i] = exp(-r*T)*(ST-K);
           		P[i] = 0;
        	}
        	ave[i] += C[i];														//Accumulate the values
			avep[i] += P[i];       
    	}
		ave[i] = ave[i]/L;
		ave2[i] = pow(ave[i],2);
		
		avep[i] = avep[i]/L;
		ave2p[i] = pow(avep[i],2);
	}
	
	output_final.open("final.dat");
	for (int i = 0; i < steps; i++) {    		//Blocking method
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
			
			sump_prog[i] += avep[j];
			sum2p_prog[i] += ave2p[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		
		sump_prog[i] = sump_prog[i]/(i+1);
		sum2p_prog[i] = sum2p_prog[i]/(i+1);
		errp_prog[i] = error(sump_prog, sum2p_prog, i);
		output_final << sum_prog [i] << "\t" << err_prog[i] << "\t" << sump_prog [i] << "\t" << errp_prog[i] << endl; 
	}
	output_final.close();
	
	//Compute the option prices through the Brownian Motion
	for (int i = 0; i < steps; i++) {			//Initialize all the values of the price and variables for the blocking method
		S[i] = 0;
		
		sum_prog[i] = 0;						//Call-option
		sum2_prog[i] = 0;
		err_prog[i] = 0;
		ave[i] = 0;
		ave2[i] = 0;
		
		sump_prog[i] = 0;						//Put-option
		sum2p_prog[i] = 0;
		errp_prog[i] = 0;
		avep[i] = 0;
		ave2p[i] = 0;
	}
	
	for (int i = 1; i < steps; i++) {			//Generate the time steps
		t[i] = 1./(double)steps + t[i-1]; 
	}
	
	for (int c = 0; c < steps; c++) { 			//Cycle over the different blocks
		for (int i = 0; i < L; i++) {			//Cycle inside a block
			for (int k = 0; k < L; k++) {
			C[k] = 0;								//Call-option
			P[k] = 0;								//Put-option
			}
			for (int k = 1; k < steps; k++) {	
				S[k] = 0;							//Set the prices to zero every time the simulation starts again
			}
			S[0] = S0;								//Set the initial price
			for (int j = 1; j < steps; j++) { 														//Simulate the Brownian Motion
				z = rnd.Gauss(0,1);																	//Normal random variable
				S[j] = S[j-1]*exp((r-0.5*pow(sigma,2))*(t[j]-t[j-1])+sigma*z*sqrt((t[j]-t[j-1])));	//Make the GBM move
			}
			ST = S[steps-1];																		//Set the final value of S
        	if (0 > (ST-K)) {																		//Choose whether to evaluate the call or put option
        	    C[i] = 0;
        		P[i] = exp(-r*T)*(K-ST);
       		} else {
            	C[i] = exp(-r*T)*(ST-K);
           		P[i] = 0;
        	}
        	ave[c] += C[i];						//Accumulate the values
			avep[c] += P[i];       
    	}
		ave[c] = ave[c]/L;
		ave2[c] = pow(ave[c],2);
		
		avep[c] = avep[c]/L;
		ave2p[c] = pow(avep[c],2);
	}
	
	output_path.open("path.dat");
	for (int i = 0; i < steps; i++) {    		//Blocking method for computing the mean value and the variance with respective uncertainties
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
			
			sump_prog[i] += avep[j];
			sum2p_prog[i] += ave2p[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		
		sump_prog[i] = sump_prog[i]/(i+1);
		sum2p_prog[i] = sum2p_prog[i]/(i+1);
		errp_prog[i] = error(sump_prog, sum2p_prog, i);
		
		output_path << sum_prog [i] << "\t" << err_prog[i] << "\t" << sump_prog [i] << "\t" << errp_prog[i] << endl;
	}
	output_path.close();
   	rnd.SaveSeed(); 
return 0;
}


void Random_Initializer() {
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
   	return;
}

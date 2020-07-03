#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "functions.h" 

using namespace std;


int main (int argc, char *argv[]){
	
	double lambda = 1.0, mu = 0, gamma = 1.0;			//Parameters for the exp and Lorentz distributions
	ofstream output_std, output_exp, output_lorentz;	//Output ofstreams
	
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
	
	//Standard dice
	for (int k = 0; k < M; k++) {					//Initialize all the arrays to 0
		S1[k] = 0;
		S2[k] = 0;
		S10[k] = 0;
		S100[k] = 0;
	}
	
	output_std.open("output_standard.dat");
	for (int j = 0; j < M; j++) {
		N = 1;								//Set N = 1
		for (int i = 0; i < N; i++) {		//Compute S1
			x = ceil(rnd.Rannyu(0,6));		//Generate a random number between 1 and 6 (dice)
			S1[j] = S1[j] + x;
			
		}
		N = 2;								//Set N = 2
		for (int i = 0; i < N; i++) {		//Compute S2
			x = ceil(rnd.Rannyu(0,6));		//Generate a random number between 1 and 6 (dice)
			S2[j] = S2[j] + x;
		}
		N = 10;								//Set N = 10
		for (int i = 0; i < N; i++) {		//Compute S10
			x = ceil(rnd.Rannyu(0,6));		//Generate a random number between 1 and 6 (dice)
			S10[j] = S10[j] + x;
		}
		N = 100;							//Set N = 100
		for (int i = 0; i < N+1; i++) {		//Compute S100
			x = ceil(rnd.Rannyu(0,6));		//Generate a random number between 1 and 6 (dice)
			S100[j] = S100[j] + x;
		}
		output_std << S1[j] << "\t" << S2[j]/2 << "\t"<< S10[j]/10 << "\t"<< S100[j]/100 << endl;
	}
	output_std.close();
	
	//Exponential distribution
	for (int k = 0; k < M; k++) {		//Reset the arrays for the S's
		S1[k] = 0;
		S2[k] = 0;
		S10[k] = 0;
		S100[k] = 0;
	}
	
	output_exp.open("output_exp.dat");			//Same as for the standard dice
	for (int j = 0; j < M; j++) {
		N = 1;
		for (int i = 0; i < N; i++) {
			x = (rnd.Exponential(lambda));		//Generate an exponential random number
			S1[j] = S1[j] + x;
		}
		N = 2;
		for (int i = 0; i < N; i++) {
			x = (rnd.Exponential(lambda));
			S2[j] = S2[j] + x;
		}
		N = 10;
		for (int i = 0; i < N; i++) {
			x = (rnd.Exponential(lambda));
			S10[j] = S10[j] + x;
		}
		N = 100;
		for (int i = 0; i < N; i++) {
			x = (rnd.Exponential(lambda));
			S100[j] = S100[j] + x;
		}
		output_exp << S1[j] << "\t" << S2[j]/2 << "\t"<< S10[j]/10 << "\t"<< S100[j]/100 << endl;
	}
	output_exp.close();
	
	//Exponential distribution. Same as the distributions before
	for (int k = 0; k < M; k++) {
		S1[k] = 0;
		S2[k] = 0;
		S10[k] = 0;
		S100[k] = 0;
	}

	output_lorentz.open("output_lorentz.dat");
	for (int j = 0; j < M; j++) {
		N = 1;
		for (int i = 0; i < N; i++) {
			x = (rnd.Lorentz(mu,gamma)); 		//Generate a Lorentzian random number 
			S1[j] = S1[j] + x;
		}
		N = 2;
		for (int i = 0; i < N; i++) {
			x = (rnd.Lorentz(mu,gamma));
			S2[j] = S2[j] + x;
		}
		N = 10;
		for (int i = 0; i < N; i++) {
			x = (rnd.Lorentz(mu,gamma));
			S10[j] = S10[j] + x;
		}
		N = 100;
		for (int i = 0; i < N; i++) {
			x = (rnd.Lorentz(mu,gamma));
			S100[j] = S100[j] + x;
		}
		output_lorentz << S1[j] << "\t" << S2[j]/2 << "\t"<< S10[j]/10 << "\t"<< S100[j]/100 << endl;
	}
	output_lorentz.close();
	rnd.SaveSeed();
	
return 0;
}

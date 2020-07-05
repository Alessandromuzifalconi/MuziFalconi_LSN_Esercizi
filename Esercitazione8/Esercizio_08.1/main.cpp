#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "functions.h"


using namespace std;

int main (int argc, char *argv[]){
	
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
	
	n_acc = 0.0;
	n_rej = 0.0;
	
	//Starting position 
	x_old = 2.5;
	//Starting parameters
	mu = 1;
	sigma = sqrt(0.8);
	
	ofstream x_pos("x_pos.dat");
	//Generate the positions distributed as |Psi|^2
	for (int i = 0; i < M; i++) {	
		x_new = x_old + rnd.Rannyu(-x_max, x_max);	 					//Move the particle
		//Implementing the minimum
		alpha = 1;						
		if (pow(Psi(x_new, sigma, mu),2)/pow(Psi(x_old, sigma,mu),2) < 1)
			alpha = pow(Psi(x_new, sigma, mu),2)/pow(Psi(x_old, sigma, mu),2);	
		if (rnd.Rannyu() <= alpha) {
			x_old = x_new;
			n_acc = n_acc+1;
		} else {
			n_rej = n_rej+1;
		}
		x[i] =  x_old;	
		x_pos << x[i] << endl;
	}		
	x_pos.close();
	cout << "Accettati/rigettati \t" << n_acc/n_rej << endl;
		
	//Compute expectation value <H> with blocking method
	for (int i = 0; i < N; i++) {								//Cycle over blocks
		sum = 0;
		for (int j = 0; j < L; j++) {							//Cycle inside a block
			k = j + i*L;
			sum += 	HPsi(x[k],sigma,mu)/Psi(x[k], sigma, mu);
		}
		ave[i] = sum/L;
		ave2[i] = pow(ave[i],2);
	}
	
	ofstream H("H.dat");
	for (int i = 0; i < N; i++) {    		//Blocking method for computing the mean value and the variance with respective uncertainties
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		err_prog[i] = Error(sum_prog, sum2_prog, i);
		H << sum_prog [i] << "\t" << err_prog[i] << endl; 			//Write results on a file
	}
	H.close();
return 0;
}

double Psi (double x, double sigma, double mu) {
	return exp(-pow(x-mu,2)/(2*pow(sigma,2))) + exp(-pow(x+mu,2)/(2*pow(sigma,2)));
}

double V (double x) {
	return pow(x,4)-2.5*pow(x,2);
}

double Error (double* ave, double* ave2, int i) {
	 if( i == 1 ) return 0.0;
	return sqrt((ave2[i]-pow(ave[i],2))/i);
}

double HPsi(double x, double sigma, double mu) {
	double e1 = exp(-pow(x - mu, 2) / (2 * pow(sigma,2)));
	double e2 = exp(-pow(x + mu, 2) / (2 * pow(sigma,2)));
	return  0.5*1/pow(sigma,2)*e1 - pow(x-mu,2)/pow(sigma,4) *e1 + 1/pow(sigma,2)*e2 - pow(x + mu, 2) / pow(sigma,4)*e2 + V(x) * Psi(x, sigma, mu);
}

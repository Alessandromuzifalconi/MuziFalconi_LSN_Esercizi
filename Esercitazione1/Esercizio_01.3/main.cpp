#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistics.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

	
	ofstream results;
	
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
				
	for (int i = 0; i < N; i++) { 			//Cycle over different blocks
		hits = 0;
		for (int j = 0; j < K; j++) {		//Cycle inside a block
			x = rnd.Rannyu(0,2*d);			//Throw the needle between to lines. x is the x-coordinate of the left tip of the needle
			theta = rnd.Theta_rejection();	//Generate a random angle in [0,Pi/2] using the rejection techinque
			l = L*cos(theta);				//Project the needle along the x-axis
			if (x < d && x+l > d) {			//Check if the needle crosses the first line
				hits = hits+1;
			}
			if (x > d && x + l > 2*d) {		//Check if the needle crosses the second line
				hits = hits+1;
			}
		}
		pi[i] = K*2*L/double(d*hits);
		pi2[i] = pow(pi[i],2);
	}	
	
	results.open("Buffon_output.dat");	
	for (int i = 0; i < N; i++) {    		//Blocking method for computing pi and uncertainty
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += pi[j];
			sum2_prog[i] += pi2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		results << sum_prog[i] << "\t" << err_prog[i] << endl;
	}
	results.close();
	rnd.SaveSeed();
return 0;
}

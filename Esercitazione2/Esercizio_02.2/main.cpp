#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistics.h"
#include "position.h"
#include "functions.h"

using namespace std;


int main (int argc, char *argv[]){

	ofstream output_discr, output_cont;
	
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
	
	//Discrete random walk	
	for (int i = 0; i < N; i++) { 				//Cycle over the different repetitions of the random walk
		pos.SetX(0);							//Set the starting position of the RW in (0,0,0)
		pos.SetY(0);
		pos.SetZ(0);
		for (int k = 0; k < steps; k++) {	
				r[k] = 0;						//Set the distances to zero every time the rw starts again
		}
		for (int j = 0; j < steps; j++) {
			x = rnd.Rannyu(); 				
			if (x >= 0.5) {						//If x>=0.5 make a step forward 
				x = rnd.Rannyu();
					if (x < 0.33) 				
						pos.SetX(pos.GetX()+a);	//Move in the x direction
					if (x < 0.66) 
						pos.SetY(pos.GetY()+a); //Move in the y direction
					else
						pos.SetZ(pos.GetZ()+a); //Move in the z direction
			}
			else { 								//If x<0.5 make a step backwards
				x = rnd.Rannyu();
					if (x < 0.33) 				//Move in the x direction
						pos.SetX(pos.GetX()-a);
					if (x < 0.66) 
						pos.SetY(pos.GetY()-a); //Move in the y direction
					else
						pos.SetZ(pos.GetZ()-a); //Move in the z direction
			}
			r[j] = r[j] + pos.GetDist();		//Distance from the origin at step j
			sums[j] += pow(r[j],2);				//Accumulate the distances at step j
		}
	}
	
	for (int j = 0; j < steps; j++) {			
		sums[j] = sums[j]/N;					//Divide by the number of rw made
		sums[j] = sqrt(sums[j]);
		sums2[j] = pow(sums[j],2);				
	}
	output_discr.open("RW_discr");
	for (int i = 0; i < steps; i++) {    		//Blocking method 
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += sums[j];
			sum2_prog[i] += sums2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		output_discr << sum_prog [i] << "\t" << err_prog[i] << endl; 		
	}
	output_discr.close();
	
	//Continuos random walk
	for (int i = 0; i < steps; i++) {			//Iinitialize all the positions and variables for the blocking method to 0
		r[i] = 0;
		sum_prog[i] = 0;
		sum2_prog[i] = 0;
		err_prog[i] = 0;
		sums[i] = 0;
		sums2[i] = 0;
	}
	
   	for (int i = 0; i < N; i++) { 				//Cycle over the different repetitions of the random walk
		pos.SetX(0);							//Set the starting position of the RW in (0,0,0)
		pos.SetY(0);
		pos.SetZ(0);
		for (int k = 0; k < steps; k++) {		//Set the distances to zero every time the rw starts again
				r[k] = 0;
		}
		for (int j = 0; j < steps; j++) {
			theta = rnd.Rannyu(0,M_PI);						//Generate the random angles for the direction
			phi = rnd.Rannyu(0,2*M_PI);
			pos.SetX(pos.GetX()-a*sin(theta)*cos(phi));		//Move 
			pos.SetY(pos.GetY()-a*sin(theta)*sin(phi));
			pos.SetZ(pos.GetZ()-a*cos(theta));				
			r[j] = r[j] + pos.GetDist();					//Distance from the origin at step j
			sums[j] += pow(r[j],2);							//Accumulate the distances at step j
		}
	}
	
	for (int j = 0; j < steps; j++) {
		sums[j] = sums[j]/N;
		sums[j] = sqrt(sums[j]);
		sums2[j] = pow(sums[j],2);
	}
	
	output_cont.open("RW_cont");
	for (int i = 0; i < steps; i++) {    		//Blocking method
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += sums[j];
			sum2_prog[i] += sums2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		output_cont << sum_prog [i] << "\t" << err_prog[i] << endl; 	
	}
	output_cont.close();
   	rnd.SaveSeed(); 
return 0;
}

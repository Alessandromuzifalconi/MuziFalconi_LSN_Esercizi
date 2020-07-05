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
	
	for(int k = 0; k < steps; k++){
		ave[k] = new double[nblocks] {0};
		ave2[k] = new double[nblocks] {0};
  	}
  	
	
	//Discrete random walk	
	cout << "Discrete RW" << endl;
	for (int i = 0; i < nblocks; i++) { 		//Cycle over different blocks 
		for (int k = 0; k < L; k++) {			//Cycle inside a block
			pos.SetX(0);						//Set the starting position of the RW in (0,0,0)
			pos.SetY(0);
			pos.SetZ(0);
			for (int j = 1; j < steps; j++) {
				x = floor(rnd.Rannyu(0,2)); 
				if (x == 1) {						//If x==1 make a step forward 
					x = floor(rnd.Rannyu(0,3));
						if (x == 0) 				
							pos.SetX(pos.GetX()+a);	//Move in the x direction
						else if (x == 1) 
							pos.SetY(pos.GetY()+a); //Move in the y direction
						else
							pos.SetZ(pos.GetZ()+a); //Move in the z direction
				} else { 								//If x<0.5 make a step backwards
					x = floor(rnd.Rannyu(0,3));
						if (x == 0) 				
							pos.SetX(pos.GetX()-a);
						else if (x == 1) 
							pos.SetY(pos.GetY()-a); //Move in the y direction
						else
							pos.SetZ(pos.GetZ()-a); //Move in the z direction
				}
				ave[j][i] += pow(pos.GetDist(),2)/double(L);	//Accumulate the distances: j = step, i = block
				ave2[j][i] = pow(ave[j][i],2);
			}
		}	
	}
		
	output_discr.open("RW_discr");
	for (int j = 0; j < steps; j++) {    		//Blocking method 
		sum_prog[j] = 0;
		sum2_prog[j] = 0;
		for (int i = 0; i < nblocks; i++) {
			sum_prog[j] += ave[j][i];
			sum2_prog[j] += ave2[j][i];
		}	
		sum_prog[j]/=nblocks;
		sum2_prog[j]/=nblocks;
		err_prog[j] = error(sum_prog[j], sum2_prog[j], nblocks);
		output_discr << sum_prog[j] << "\t" << err_prog[j]  << endl; 
	}
			
	output_discr.close();
	
	//Continuos random walk
	cout << "Continuos RW" << endl;
	for (int j = 0; j < steps; j++) {			//Iinitialize all the positions and variables for the blocking method to 0
		sum_prog[j] = 0;
		sum2_prog[j] = 0;
		err_prog[j] = 0;
		for (int i = 0; i < nblocks; i++) {
			ave[j][i] = 0;	
			ave2[j][i] = 0;
		}			
	}
	
   for (int i = 0; i < nblocks; i++) { 		//Cycle over different blocks 
		for (int k = 0; k < L; k++) {			//Cycle inside a block
			pos.SetX(0);							//Set the starting position of the RW in (0,0,0)
			pos.SetY(0);
			pos.SetZ(0);
			for (int j = 1; j < steps; j++) {
				theta = rnd.Rannyu(0,M_PI);						//Generate the random angles for the direction
				phi = rnd.Rannyu(0,2*M_PI);
				pos.SetX(pos.GetX()-a*sin(theta)*cos(phi));		//Move 
				pos.SetY(pos.GetY()-a*sin(theta)*sin(phi));
				pos.SetZ(pos.GetZ()-a*cos(theta));				
				
				ave[j][i] += pow(pos.GetDist(),2)/double(L);	//Accumulate the distances: j = step, i = block
				ave2[j][i] = pow(ave[j][i],2);
			}
		}
	}
	
	output_cont.open("RW_cont");
	for (int j = 0; j < steps; j++) {    		//Blocking method 
		sum_prog[j] = 0;
		sum2_prog[j] = 0;
		for (int i = 0; i < nblocks; i++) {
			sum_prog[j] += ave[j][i];
			sum2_prog[j] += ave2[j][i];
		}	
		sum_prog[j]/=nblocks;
		sum2_prog[j]/=nblocks;
		err_prog[j] = error(sum_prog[j], sum2_prog[j], nblocks);
		output_cont << sum_prog [j] << "\t" << err_prog[j] << endl; 	
	}
	output_cont.close();
	
   	rnd.SaveSeed(); 
return 0;
}

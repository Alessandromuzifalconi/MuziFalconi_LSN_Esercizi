#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "position.h"
#include "statistics.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){
	
	ofstream psi100_unif, psi100_gauss, psi210_unif, psi210_gauss;
	ofstream psi100_u_conf, psi100_g_conf, psi210_u_conf, psi210_g_conf;
	
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
   	
	// Psi100 uniform
	r_max = 1.2;												//So that we have a 50% acceptance rate
	
	x1.SetX(5);													//Set starting position
	x1.SetY(5);
	x1.SetZ(5);
	
	psi100_u_conf.open("Psi100_u_conf.dat");
	for (int i = 0; i < N; i++) {								//Cycle over blocks
		for (int j = 0; j < L; j++) {							//Cycle inside a block
			x2.SetX(x1.GetX() + rnd.Rannyu(-r_max, r_max));		//Move at random around the point (uniform distribution)
			x2.SetY(x1.GetY() + rnd.Rannyu(-r_max, r_max));
			x2.SetZ(x1.GetZ() + rnd.Rannyu(-r_max, r_max));

			//Implementing the minimum
			alpha = 1;							
			if (Psi100(x2)/Psi100(x1) < 1) 			//If probability density of x2 < prob. dens. x1 set alpha < 1
				alpha = Psi100(x2)/Psi100(x1);
		
			if (rnd.Rannyu() <= alpha) { 		 	//Decide whether to accept or reject the move		
				x1.SetX(x2.GetX());
				x1.SetY(x2.GetY());
				x1.SetZ(x2.GetZ());
				n_acc = n_acc + 1;				
			} else {	
				n_rej = n_rej+1;
			}
			if (j%nsave == 0)							//Save data for plotting the points 	
				psi100_u_conf << x1.GetX() << "\t" << x1.GetY() << "\t" << x1.GetZ() << endl;
			ave[i] +=  x1.GetDist();				//Accumulate the distances
		}	
	}	
	psi100_u_conf.close();
	
	for (int j = 0; j < N; j++) {					//Divide for the number of elements in a block
		ave[j] = ave[j]/L;					
		ave2[j] = pow(ave[j],2);				
	}	
	
	psi100_unif.open("Psi100_u.dat");
	for (int i = 0; i < N; i++) {    		//Blocking method for computing the mean value and the variance with respective uncertainties
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		psi100_unif << sum_prog [i] << "\t" << err_prog[i] << endl; 			//Write results on a file
	}
	psi100_unif.close();
	
	cout << "\t Psi100 uniform: " << endl/* << " Media: " << mean(M, r) << endl*/;
	cout << " Accepted/Rejected: " << (double)n_acc/n_rej << endl << endl;
	
	// Psi100 Gauss
	n_acc = 0;
	n_rej = 0;
	
	r_max = 0.75;
	
	x1.SetX(5);
	x1.SetY(5);
	x1.SetZ(5);
	
	for (int i = 0; i < N; i++) {
		ave[i] = 0;
		ave2[i] = 0;
		sum_prog[i] = 0;
		sum2_prog[i] = 0;
		err_prog[i] = 0;
	}
	
	psi100_g_conf.open("Psi100_g_conf.dat");
	for (int i = 0; i < N; i++) {							//Cycle over blocks
		for (int j = 0; j < L; j++) {						//Cycle inside a block
			x2.SetX(rnd.Gauss(x1.GetX(), (r_max)));		//Move around X1 with a gaussian distribution
			x2.SetY(rnd.Gauss(x1.GetY(), (r_max)));
			x2.SetZ(rnd.Gauss(x1.GetZ(), (r_max)));

			alpha = 1;
			if (Psi100(x2)/Psi100(x1) < 1) 
				alpha = Psi100(x2)/Psi100(x1);
		
			if (rnd.Rannyu() <= alpha) {
				x1.SetX(x2.GetX());
				x1.SetY(x2.GetY());
				x1.SetZ(x2.GetZ());
				n_acc = n_acc+1;
			} else {
				n_rej = n_rej+1;
			}
			if (j%nsave == 0)							//Save data for plotting the points 	
				psi100_g_conf << x1.GetX() << "\t" << x1.GetY() << "\t" << x1.GetZ() << endl;
			ave[i] +=  x1.GetDist();
		}	
	}
	psi100_g_conf.close();
	
	for (int j = 0; j < N; j++) {			
		ave[j] = ave[j]/L;					
		ave2[j] = pow(ave[j],2);				
	}	
	
	psi100_gauss.open("Psi100_G.dat");
	for (int i = 0; i < N; i++) {    		//blocking method for computing the mean value and the variance with respective uncertainties
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		psi100_gauss << sum_prog [i] << "\t" << err_prog[i] << endl; 			//write results on a file
	}
	psi100_gauss.close();
	
	cout << "\t Psi100 gauss: " << endl/* << " Media: " << mean(M, r) << endl*/;
	cout << " Accepted/Rejected: " << (double)n_acc/n_rej << endl << endl;
	
	// Psi210 unif
	n_acc = 0;
	n_rej = 0;
	
	r_max = 2.9;
	
	x1.SetX(8);
	x1.SetY(8);
	x1.SetZ(8);
	
	
	for (int i = 0; i < N; i++) {
		ave[i] = 0;
		ave2[i] = 0;
		sum_prog[i] = 0;
		sum2_prog[i] = 0;
		err_prog[i] = 0;
	}
	psi210_u_conf.open("Psi210_u_conf.dat");
	for (int i = 0; i < N; i++) {								//Cycle over blocks
		for (int j = 0; j < L; j++) {							//Cycle inside a block
			x2.SetX(x1.GetX() + rnd.Rannyu(-r_max, r_max));		//Move at random around the point (uniform distribution)
			x2.SetY(x1.GetY() + rnd.Rannyu(-r_max, r_max));
			x2.SetZ(x1.GetZ() + rnd.Rannyu(-r_max, r_max));

			//Implementing the minimum
			alpha = 1;							
			if (Psi210(x2)/Psi210(x1) < 1) 			//If probability density of x2 < prob. dens. x1 set alpha < 1
				alpha = Psi210(x2)/Psi210(x1);
		
			if (rnd.Rannyu() <= alpha) { 		 	//Decide whether to accept or reject the move		
				x1.SetX(x2.GetX());
				x1.SetY(x2.GetY());
				x1.SetZ(x2.GetZ());
				n_acc = n_acc + 1;				
			} else {	
				n_rej = n_rej+1;
			}
			if (j%nsave == 0)							//Save data for plotting the points 	
				psi210_u_conf << x1.GetX() << "\t" << x1.GetY() << "\t" << x1.GetZ() << endl;
			ave[i] +=  x1.GetDist();				//Accumulate the distances
		}	
	}
	psi210_u_conf.close();
	
	for (int j = 0; j < N; j++) {					//Divide for the number of elements in a block
		ave[j] = ave[j]/L;					
		ave2[j] = pow(ave[j],2);				
	}	
	
	psi210_unif.open("Psi210_u.dat");
	for (int i = 0; i < N; i++) {    		//Blocking method for computing the mean value and the variance with respective uncertainties
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		psi210_unif << sum_prog [i] << "\t" << err_prog[i] << endl; 			//Write results on a file
	}
	psi210_unif.close();
	
	cout << "\t Psi210 uniform: " << endl;
	cout << " Accepted/Rejected: " << (double)n_acc/n_rej << endl << endl;
	
	// Psi210 Gauss
	n_acc = 0;
	n_rej = 0;
	
	r_max = 1.9	;
	
	// Starting position is (4,4,4) in units of Bohr radius
	x1.SetX(8);
	x1.SetY(8);
	x1.SetZ(8);
	
	for (int i = 0; i < N; i++) {
		ave[i] = 0;
		ave2[i] = 0;
		sum_prog[i] = 0;
		sum2_prog[i] = 0;
		err_prog[i] = 0;
	}
	
	psi210_g_conf.open("Psi210_g_conf.dat");
	for (int i = 0; i < N; i++) {								//Cycle over blocks
		for (int j = 0; j < L; j++) {							//Cycle inside a block
			x2.SetX(rnd.Gauss(x1.GetX(), (r_max)));				//Move around X1 with a gaussian distribution
			x2.SetY(rnd.Gauss(x1.GetY(), (r_max)));
			x2.SetZ(rnd.Gauss(x1.GetZ(), (r_max)));

			alpha = 1;
			if (Psi210(x2)/Psi210(x1) < 1) 
				alpha = Psi210(x2)/Psi210(x1);
		
			if (rnd.Rannyu() <= alpha) {
				x1.SetX(x2.GetX());
				x1.SetY(x2.GetY());
				x1.SetZ(x2.GetZ());
				n_acc = n_acc+1;
			} else {
				n_rej = n_rej+1;
			}
			if (j%nsave == 0)							//Save data for plotting the points 	
				psi210_g_conf << x1.GetX() << "\t" << x1.GetY() << "\t" << x1.GetZ() << endl;
			ave[i] +=  x1.GetDist();
		}	
	}
	psi210_g_conf.close();
	
	for (int j = 0; j < N; j++) {			
		ave[j] = ave[j]/L;					
		ave2[j] = pow(ave[j],2);				
	}	
	
	psi210_gauss.open("Psi210_G.dat");
	for (int i = 0; i < N; i++) {    		//blocking method for computing the mean value and the variance with respective uncertainties
		for (int j = 0; j < i+1; j++) {
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1);
		sum2_prog[i] = sum2_prog[i]/(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
		psi210_gauss << sum_prog [i] << "\t" << err_prog[i] << endl; 			//write results on a file
	}
	psi210_gauss.close();
	
	cout << "\t Psi210 gauss: " << endl;
	cout << " Accepted/Rejected: " << (double)n_acc/n_rej << endl << endl;	
	
	return 0;
}

double Psi100 (Position r) {
	return exp(- 2*r.GetDist())/M_PI;
}

double Psi210 (Position r) {
	return 1./64*2/(M_PI)*pow(r.GetZ(),2)*exp(-r.GetDist());
}




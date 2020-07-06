#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "statistics.h"
#include "individual.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

	Individual I_old;	
	Individual I_new;			
	ofstream cities_coord, best_individual, Hall_of_Fame;
	RandomInitializer();				//Initialize the pseudo random numbers generator

	cout << endl << "Simulated annealing for solving the travelling salesman problem." << endl << "-----------------------------------------------------------------" << endl;
	//Create the first generation
	Cities_generator(method, N); 			//Create the position of the cities
   	cities_coord.open("Cities.dat");		//Write cities' positions on file
   	for (int i = 0; i < N; i++) {
   		cities_coord << cities[i][0] << "\t" << cities[i][1] << endl;
    }
    cities_coord.close();
    
	while (i_s < N) {				//Generate a random sequence
		seq_old[i_s] = (int) floor(rnd.Rannyu(1,32));
	 	j_s = 1;
	   	while (j_s<i_s) { 										//Check if the city is already in the sequence
	   		if (seq_old[j_s] == seq_old[i_s]) { 				//If the city is already in the sequence generate the city again
	   			j_s = i_s; 										//Exit from the while (j_s < i_s) in order to generate the city again
	   		   	i_s = i_s-1;
	   			if (i_s<0) 
	   				cout << "Error: i_s < 0" << endl;
	   		} else { 
	   			j_s = j_s+1;
	   		}	 
	   	}
	   	i_s = i_s+1;
	}
    I_old.SetSequence(seq_old);									//Set the order of the cities
 
   	Hall_of_Fame.open("Hall_of_Fame.dat");
   	for (int i = 0; i < m; i++) {				//Cycle over different temperatures 
		beta += (1./Tf)/double(m);				//Decrease temperature
		for (int j = 0; j < n; j++) {			//Make n MC steps for each temperature
			I_new.SetSequence(I_old.GetSequence());	
  			Exchange (I_new.GetSequence());				//Try to swap two cities
			//Implementing the minimum for Metropolis
			alpha = 1;							
			if (L(I_old.GetSequence())/L(I_new.GetSequence()) < 1) 		//If the new order is worse than the old
				alpha = exp(-beta*(L(I_new.GetSequence())-L(I_old.GetSequence())));
		
			if (rnd.Rannyu() <= alpha) { 		 	//Decide whether to accept or reject the move		
				I_old.SetSequence(I_new.GetSequence());
				n_acc = n_acc + 1;				
			} else {	
				n_rej = n_rej+1;
			}
		}	 
		Hall_of_Fame << i << "\t" << beta << "\t" << L(I_old.GetSequence()) << endl; //Write on file what has happened in this generation
	}
	Hall_of_Fame.close();
    best_individual.open("Best.dat");		//Write on file the sequence of cities at the end of the annealing 
    for (int i = 0; i < N; i++) { 
    	best_individual << I_old.GetGene(i) << "\t" << cities[I_old.GetGene(i)][0] << " " << cities[I_old.GetGene(i)][1] << endl;
    }
    best_individual.close();
   	rnd.SaveSeed(); 	

return 0;
}


void RandomInitializer() {
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
}
	
void Cities_generator (int method, int N) {
	if (method == 0) { 						//Generate the cities on a circumference (R=1)
		cout << "Travel in " << N << " cities placed randomly on a circumference of radius 1." << endl;
		double theta = 0;
		for (int i = 0; i < N; i++) {
			theta = rnd.Theta();
			cities[i][0] = cos(theta); 	//x
			cities[i][1] = sin(theta); 	//y
		}
	} else { 
		cout << "Travel in " << N << " cities placed randomly inside a square of size 2." << endl;
		for (int i = 0; i < N; i++) {
			cities[i][0] = rnd.Rannyu(-1.,1.); 	//x
			cities[i][1] = rnd.Rannyu(-1.,1.); 	//y
		}
	}	
}	

int PBC(int i) { 						//Periodic Boundary Condition
        return i - 32 * rint(i/32); 	//rint rounds x to an integral value
}

int Check(int* sequence) { 				//Check if a city appears twice in a sequence. If there is a problem returns 1
    int control = 0;
    for (int i = 0; i < N; i++) {			
        for (int j = i+1; j < N; j++) {
            if (sequence[j] == sequence[i]) {
                control = 1;
            }
        }
    }
    return control;
}

double L (int* sequence) {
	double sum = 0;
	for (int i = 0; i < N; i++) {
		sum += sqrt(pow(cities[sequence[PBC(i)]][0]-cities[sequence[PBC(i+1)]][0],2) + pow(cities[sequence[PBC(i)]][1]-cities[sequence[PBC(i+1)]][1],2));

	}
	return sum;
}

void Exchange (int* sequence) {
	int* appo = new int [N];
	
	for (int i = 0; i < N; i++) 
		appo[i] = sequence[i];
	
	int pos1 = (int) floor(rnd.Rannyu(1.,32.));
	int pos2 = (int) floor(rnd.Rannyu(1.,32.));    	//Randomly select two positions to be swapped
	int gene1 = sequence[pos1];
	int gene2 = sequence[pos2];
	sequence[pos1] = gene2;
	sequence[pos2] = gene1;
	if (Check(sequence) == 1) {						 //If there is a problem undo the swap
		for (int i = 0; i < N; i++) 
			sequence[i] = appo[i];
	}
	return;
}

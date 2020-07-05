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

	Individual I[P];					//Population
	ofstream cities_coord, best_individual, Hall_of_Fame;
	RandomInitializer();				//Initialize the pseudo random numbers generator

	cout << endl << "Genetic algorithm for solving the travelling salesman problem." << endl << "-----------------------------------------------------------------" << endl;
	//Create the first generation
	Cities_generator(method, N); 			//Create the position of the cities
   	cities_coord.open("Cities.dat");		//Write cities' positions on file
   	for (int i = 0; i < N; i++) {
   		cities_coord << cities[i][0] << "\t" << cities[i][1] << endl;
    }
    cities_coord.close();
    cout << "Create the first generation of a population of " << P << " individuals." 	<< endl << "----------------------------------------------------------------" << endl;
    Population_generator(I);				//Create the population
	Order_Individuals(I,order);				//Order individuals form best to worst
	
	
   	Hall_of_Fame.open("Hall_of_Fame.dat");
   	for (int g = 1; g <= Generations; g++) {	//Cycle over generations 
		
		Mutate(I);
		Order_Individuals(I,order);		//Order individuals form best to worst
		
		for (int k = 0; k < P; k++) {	//Save distances of each individual
			lenghts[k] = L(I[k].GetSequence());
		} 
		
		if (g%(Generations/100) == 0)	//Print on screen
			cout << "Generation: " << g << ".\t The best individual was in position: " << order[0] <<".\t And has lenght: " << L(I[0].GetSequence()) << endl;
		
		Hall_of_Fame << g << "\t" << L(I[0].GetSequence()) << "\t" << mean(P/2,lenghts) << "\t" << mean(P,lenghts) << endl; //Write on file what has happened in this generation
	}
	Hall_of_Fame.close();
   
    best_individual.open("Best.dat");		//Write on file the charachteristics of the best individual of the last generation
    for (int i = 0; i < N; i++) { 
    	best_individual << I[0].GetGene(i) << "\t" << cities[I[0].GetGene(i)][0] << " " << cities[I[0].GetGene(i)][1] << endl;
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

void Population_generator (Individual* I) {
	int random_seq[N];
	int i_s = 1, j_s = 1;
	for (int n = 0; n < P; n++) { 			//Cycle over different indivdiuals inside the population
	   		i_s = 1;
	   		random_seq[0] = 0; 				//The first city is fixed at 0
	   		while (i_s < N) {				//Generate a random sequence
	   			random_seq[i_s] = (int) floor(rnd.Rannyu(1,32));
	   			j_s = 1;
	   			while (j_s<i_s) { 										//Check if the city is already in the sequence
	   				if (random_seq[j_s] == random_seq[i_s]) { 			//If the city is already in the sequence generate the city again
	   				   	j_s = i_s; 										//Exit from the while (j_s<i_s) in order to generate the city again
	   				   	i_s = i_s-1;
	   					if (i_s<0) 
	   						cout << "Error: i_s < 0" << endl;
	   				} else { 
	   					j_s = j_s+1;
	   				}	 
	   			}
	   			i_s = i_s+1;
	   		}
	   		I[n].SetSequence(random_seq);
	   	}
}

int PBC(int i) { 						//Periodic Boundary Condition
        return i - 31 * rint(i/32); 	//rint rounds x to an integral value
}

int Check(Individual* I, int k) { 			//Check if a city appears twice in an individual. If there is a problem returns 1. To be used inside an individual, not over different individuals
    int control = 0;
    for (int i = 0; i < N; i++) {			//Cycle over different genes (i.e. cities)	
        for (int j = i+1; j < N; j++) {
            if (I[k].GetGene(j) == I[k].GetGene(i)) {
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

void L_Calculator (Individual* I, double* lenghts) {
	for (int i = 0; i < P; i++) 
		lenghts[i] = L(I[i].GetSequence());
	return;
}

void Order_Individuals (Individual* I, int* order) {
	Individual MrWolf[P];			
   	double* l = new double [P] {0};	
   	L_Calculator(I,l); 
    for (int i = 0; i < P; i++) {  //At the beginning the order of the individuals is just 0,...,P-1
		order[i] = i;
	}
    Quicksort(l,order,0,P-1);		//Order the distances
   	for (int i = 0; i < P; i++) {
   		MrWolf[i].SetSequence(I[order[i]].GetSequence());
	}
	for (int i = 0; i < P; i++) {
   		I[i].SetSequence(MrWolf[i].GetSequence());
	}
}

int Selection (int j) {
	if (rnd.Rannyu() < 0.5*(pow(abs(P-j+1),3))/(double)pow(P,3))
		return 1;
	else 
		return 0;
}

int Selection_Mutation (int j) {
	if (rnd.Rannyu() < 0.3*(pow(abs(P-j+1),3))/(double)pow(P,3))
		return 1;
	else 
		return 0;
}

void Mutation_PairExchange (Individual* I, int j) {
	Individual MrWolf;
	
	MrWolf.SetSequence(I[j].GetSequence());
	int pos1 = (int) floor(rnd.Rannyu(1.,32.));
	int pos2 = (int) floor(rnd.Rannyu(1.,32.));    //Randomly select two positions to be swapped
	int gene1 = I[j].GetGene(pos1);
	int gene2 = I[j].GetGene(pos2);	
	I[j].SetGene(pos1, gene2);
	I[j].SetGene(pos2, gene1);
	if (Check(I,j) == 1) {						 //If there is a problem undo the mutation
		I[j].SetSequence(MrWolf.GetSequence());
	}
	return;
}

void Mutation_Shift (Individual* I, int j) {
	Individual MrWolf;
	
	MrWolf.SetSequence(I[j].GetSequence());
    int start = (int)floor(rnd.Rannyu(1., 32.));		//Select the first gene to be shifted: possible values [1,31]
    int n = (int)ceil(rnd.Rannyu(1., 15.)); 			//Select a random n = size of the shift
    int m = (int)ceil(rnd.Rannyu(1., 15.)); 			//Select a random m = number of genes to be shifted
    int* appo_m = new int[m];
    int* appo_n = new int[n];
    for (int k = 0; k < m; k++) {						//Fill appo_m with the genes to be shifted
        appo_m[k] = I[j].GetGene(PBC(start + k)); 
    }
    for (int k = 0; k < n; k++) {
        appo_n[k] = I[j].GetGene(PBC(start + m + k)); 	//Fill appo_n with the genes to be replaced and moved
   	}
    for (int k = 0; k < n; k++) {
        I[j].SetGene(PBC(start + k), appo_n[k]);
    }
    for (int k = 0; k < m; k++) {
        I[j].SetGene(PBC(start + k + n), appo_m[k]);   
   	}
   	if (Check(I,j) == 1) {						 		//If there is a problem undo the mutation
		I[j].SetSequence(MrWolf.GetSequence());
	}
    return;
}

void Mutation_Permutation (Individual* I, int j) {
	Individual MrWolf;
	
	MrWolf.SetSequence(I[j].GetSequence());
	int start = (int)floor(rnd.Rannyu(1., 32.));		//Select the first gene: possible values [1,31]
	int m = (int)ceil(rnd.Rannyu(1., 15.)); 			//Select a random m = number of genes to be permutated	
	int* appo= new int[m];
    for (int k = start; k < m+start; k++) {				//Fill appo_m with the genes to be permutated
        appo[k-start] = I[j].GetGene(PBC(k)); 

    }
    for (int k = 0; k < m; k++) {						
        I[j].SetGene(PBC(start+k+m), appo[k]);
    }
    for (int k = 0; k < m; k++) {
        I[j].SetGene(PBC(start+k), MrWolf.GetGene(PBC(start+k+m)));
    }	
	if (Check(I,j) == 1) {						 		//If there is a problem undo the mutation
		I[j].SetSequence(MrWolf.GetSequence());
	}
    return;
}

void Mutation_Inversion (Individual* I, int j) {
	Individual MrWolf;
	
	MrWolf.SetSequence(I[j].GetSequence());
	int start = (int)floor(rnd.Rannyu(1., 32.));		//Select the first gene to be shifted: possible values [1,31]
	int m = (int)ceil(rnd.Rannyu(1., 15.)); 			//Select a random m = number of genes to be inverted	
	int* appo= new int[m];
	for (int k = start; k < m+start; k++) {				//Fill appo_m with the genes to be mutated
        appo[k-start] = I[j].GetGene(PBC(k)); 	

    }
    for (int k = 0; k < m; k++) {						
        I[j].SetGene(PBC(start+k), appo[m-k-1]);
    }
	if (Check(I,j) == 1) {						 		//If there is a problem undo the mutation
		I[j].SetSequence(MrWolf.GetSequence());
	}
    return;
}

void CrossingOver(Individual* I, int a, int b) {		//Make crossing over between to individuals

	Individual Parent1;
	Individual Parent2;
	
	Parent1.SetSequence(I[a].GetSequence());
	Parent2.SetSequence(I[b].GetSequence());
	
	int start = (int)floor(rnd.Rannyu(1., 32.));	
	int m = N-start;
	for (int i = start; i < m+start; i++) {
		I[a].SetGene(PBC(i),Parent2.GetGene(PBC(i))); 
		I[b].SetGene(PBC(i),Parent1.GetGene(PBC(i))); 
	}
	if (Check(I,a) == 1 || Check(I,b) == 1) {			//If there is a problem undo the mutation
		I[a].SetSequence(Parent1.GetSequence());
		I[b].SetSequence(Parent1.GetSequence());
	}
	return;
}

void Mutate (Individual* I) {
	for (int i = 0; i < P; i++) {			//Cycle over individuals
		if (Selection(i) == 1) {			//See if individual i is selected for crossing over
			for (int j = 0; j < P; j++) {	//Search for a partner for crossing over
				if (Selection(j) == 1 && rnd.Rannyu() < P_Crossing) {	//See if j is selected as a partner for crossing over 
					CrossingOver(I,i,j);
				}
			}
		}
		if (Selection_Mutation(i) == 1 && rnd.Rannyu() < P_PairExchange) { //See if there are mutations in individual i
			Mutation_PairExchange(I,i);
		}	
		if (Selection_Mutation(i) == 1 && rnd.Rannyu() < P_Shift) {
			Mutation_Shift(I,i);
		}
		if (Selection_Mutation(i) == 1 && rnd.Rannyu() < P_Inversion) {
			Mutation_Inversion(I,i);
		}
		if (Selection_Mutation(i) == 1 && rnd.Rannyu() < P_Inversion) {
			Mutation_Inversion(I,i);
		}
	}
	return;
}

void Quicksort(double* x, int* o, int first, int last) { 	//Function from ordering from best to worst
    int i, j, pivot;
    double temp, temp2;
    if (first < last) {
        pivot = first;
        i = first;
        j = last;
        while (i < j) {
            while (x[i] <= x[pivot] && i < last)
                i++;
            while (x[j] > x[pivot])
                j--;
            if (i < j) {
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
                temp2 = o[i];
                o[i] = o[j];
                o[j] = temp2;
            }
        }
        temp = x[pivot];
        temp2 = o[pivot];
       	x[pivot] = x[j];
        x[j] = temp;
        o[pivot] = o[j];
        o[j] = temp2;
        Quicksort(x,o,first, j - 1);
        Quicksort(x,o,j + 1, last);
    }
    return;
}

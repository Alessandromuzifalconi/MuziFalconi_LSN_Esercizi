/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() { 
	Input(); 											//Inizialization
	if (temp_array == 1) {								//Make a ladder of temperatures		
		for (int i = 0;i < 150;i++) {
	  		for(int iblk=1; iblk <= nblk; ++iblk){		//Simulation
				Reset(iblk);   							//Reset block averages
				for(int istep=1; istep <= nstep; ++istep)
				{
					Move(metro);
				  	Measure();
				  	Accumulate(); 						//Update block averages
				}
				Averages(iblk);   						//Print results for current block
			}
			ConfFinal(); 								//Write final configuration
    		temp = temp - 0.01;
    		beta = 1 / temp;
    		cout << "Now T is: "<< temp << endl;
    	}
	} else {										//Work with fixed temperature
		for(int iblk=1; iblk <= nblk; ++iblk){		//Simulation
			Reset(iblk);   							//Reset block averages
			for(int istep=1; istep <= nstep; ++istep) {
				Move(metro);
				Measure();
				Accumulate(); 						//Update block averages
			}
			Averages(iblk);   						//Print results for current block
		}
		ConfFinal(); 								//Write final configuration
	}
  	return 0;
}


void Input(void) {
	ifstream ReadInput, ReadConfig;

  	cout << "Classic 1D Ising model             " << endl;
  	cout << "Monte Carlo simulation             " << endl << endl;
  	cout << "Nearest neighbour interaction      " << endl << endl;
  	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  	cout << "The program uses k_B=1 and mu_B=1 units " << endl;


   	int p1, p2; 						//Read seed for random numbers
   	ifstream Primes("Primes");
   	Primes >> p1 >> p2 ;
   	Primes.close();

   	ifstream input("seed.in");
   	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   	rnd.SetRandom(seed,p1,p2);
   	input.close();
  

  	ReadInput.open("input.dat"); 		//Read input informations

  	ReadInput >> temp;
  	beta = 1.0/temp;
  	cout << "Temperature = " << temp << endl;

  	ReadInput >> nspin;
  	cout << "Number of spins = " << nspin << endl;

  	ReadInput >> J;
  	cout << "Exchange interaction = " << J << endl;

  	ReadInput >> h;
  	cout << "External field = " << h << endl << endl;
    
  	ReadInput >> metro; // if=1 Metropolis else Gibbs
  	ReadInput >> nblk;
  	ReadInput >> nstep;

  	if(metro==1) cout << "The program perform Metropolis moves" << endl;
  	else cout << "The program perform Gibbs moves" << endl;
  	cout << "Number of blocks = " << nblk << endl;
  	cout << "Number of steps in one block = " << nstep << endl << endl;
  
  	ReadInput >> restart;				//Read wether to restart or not
  	ReadInput >> temp_array;			//Read wether to make a ladder of temperatures between 2 and 0.5
  	ReadInput.close();


	//Prepare arrays for measurements
 	iu = 0; //Energy
  	ic = 1; //Heat capacity
  	im = 2; //Magnetization
  	ix = 3; //Magnetic susceptibility
  	n_props = 4; //Number of observables
	
	if (restart == 0) {					//Initial configuration: totally random, T-> inf
	cout << "Generate initial configuration at random" << endl;
  		for (int i=0; i<nspin; ++i) {
    		if(rnd.Rannyu() >= 0.5) s[i] = 1;
    		else s[i] = -1;
  		}
  
	//Evaluate energy etc. of the initial configuration
  	Measure();

	//Print initial values for the potential energy and virial
  	cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
	} else {						//Start from the finale configuration of the previous simulation
		cout << "Read initial configuration from file 'config.final'" << endl;
		ReadConfig.open("config.final");
		for (int i=0; i<nspin; ++i) {
			ReadConfig >> s[i];
		}
		ReadConfig.close();
	}
	return;
}
void Move(int metro) {
	int o;
	double p, energy_old, energy_new, sm, prob_up;

	for(int i=0; i<nspin; ++i) {		//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    	o = (int)(rnd.Rannyu()*nspin);	//o is the index of the particle that we are flipping
    	sm = s[o];						//sm is the value of the spin we are flipping
    	if(metro==1) { //Metropolis 
    		energy_old = Boltzmann(sm,o);
    		energy_new = Boltzmann(-sm,o);
			attempted = attempted+1;
			p = exp((-energy_new+energy_old)*beta);
			if (rnd.Rannyu() <= p) {
				s[o] = -sm;
				accepted = accepted+1;
			}
		} else { // Gibbs sampling
        	prob_up = 1./(1 + exp(-2*beta*(J*(s[Pbc(o - 1)] + s[Pbc(o + 1)]) + h)));
       		attempted = 1;
        	accepted = 1;
        	double rr = rnd.Rannyu();
        	if (rr < prob_up)
            	s[o] = +1;
        	else
            	s[o] = -1;
    	}
	}
}

double Boltzmann(int sm, int ip) { 					// Compute the energy of the new configuration
	double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  	return ene;
}

void Measure() {
  	int bin;
  	double u = 0.0, m = 0.0;

	for (int i=0; i<nspin; ++i) {											//Sum over spins
		u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);		//Energy
		m += s[i];															//Magnetization																				
  	}
  	walker[iu] = u;
  	walker[ic] = pow(u,2);
  	walker[im] = m;
  	walker[ix] = pow(m,2);
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) {				// Update block averages
	for(int i=0; i<n_props; ++i) {
    	blk_av[i] = blk_av[i] + walker[i];
   	}
   	blk_norm = blk_norm + 1.0; 		//blk_norm is counting the number of blocks done so far
}


void Averages(int iblk) //Print results for current block, iblk is the index of the block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    if (temp_array == 0) {
		cout << "Block number " << iblk << endl;
		cout << "Acceptance rate " << accepted/attempted << endl << endl;
    }
    
    if (h == 0) {
		if (metro == 1) Ene.open("output.metr.ene.0",ios::app);
		else Ene.open ("output.gibbs.ene.0",ios::app);
		stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy per spin
		glob_av[iu]  += stima_u;
		glob_av2[iu] += stima_u*stima_u;
		err_u = Error(glob_av[iu],glob_av2[iu],iblk);
		Ene << "\t" << iblk <<  "\t" << stima_u << "\t" << glob_av[iu]/(double)iblk << "\t" << err_u << endl;
		Ene.close();
    
		if (metro == 1) Heat.open("output.metr.heat.0",ios::app);
		else Heat.open ("output.gibbs.heat.0",ios::app);
		stima_c = pow(beta,2)*(blk_av[ic]/blk_norm - pow(stima_u*nspin,2))/double(nspin);
		glob_av[ic] += stima_c;
		glob_av2[ic] += stima_c*stima_c;
		err_c = Error(glob_av[ic], glob_av2[ic],iblk);
		Heat << "\t" << iblk <<  "\t" << stima_c << "\t" << glob_av[ic]/(double)iblk << "\t" << err_c << endl;
		Heat.close();
    	
    	 if (metro == 1) Chi.open("output.metr.chi.0",ios::app);
    	else Chi.open ("output.gibbs.chi.0",ios::app);
    	stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; 
    	glob_av[ix]  += stima_x;
    	glob_av2[ix] += stima_x*stima_x;
    	err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    	Chi << "\t" << iblk <<  "\t" << stima_x << "\t" << glob_av[ix]/(double)iblk << "\t" << err_x << endl;
    	Chi.close();
   	} else {
    	if (metro == 1) Mag.open("output.metr.magn.0",ios::app);
    	else Mag.open ("output.gibbs.magn.0",ios::app);
		stima_m = blk_av[im]/blk_norm/(double)nspin; 
		glob_av[im]  += stima_m;
		glob_av2[im] += stima_m*stima_m;
		err_m=Error(glob_av[im],glob_av2[im],iblk);
		Mag << "\t" << iblk <<  "\t" << stima_m << "\t" << glob_av[im]/(double)iblk << "\t" << err_m << endl;
		Mag.close();
	}
    
   
    if (temp_array == 0) 
    	cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
	if (iblk == 1) {
		return 0;
	}
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

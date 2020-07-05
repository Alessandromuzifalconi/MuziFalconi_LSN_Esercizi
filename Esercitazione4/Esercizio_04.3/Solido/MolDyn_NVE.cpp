/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
	Input();             //Inizialization
	int nconf = 1;
  	
  	int L = int(nstep/nblocks);
  	int step;
  	
  	//Prepare stuff for data blocking
	double **ave = new double*[n_props]; 
  	double **av2 = new double*[n_props];
  	
  	double sum_prog = 0, sum2_prog = 0;
  	 
  	for(int k = 0; k < n_props; k++){
		ave[k] = new double[nblocks];
		av2[k] = new double[nblocks];
  	}
	
	ofstream *ave_file = new ofstream[n_props];

  	ave_file[iv].open("ave_epot.out");
  	ave_file[ik].open("ave_ekin.out");
  	ave_file[ie].open("ave_etot.out");
  	ave_file[it].open("ave_temp.out");

	for(int iblock = 0; iblock < nblocks; ++iblock){
		
		ave[iv][iblock] = 0;		//Set averages to 0
		ave[ik][iblock] = 0;
		ave[ie][iblock] = 0;
		ave[it][iblock] = 0;
		
		for (int istep = 0; istep < L; ++istep) {	
			
			step = istep + iblock*L;
    		
    		Move();           			//Move particles with Verlet algorithm
    		
    		if(step%iprint == 0) 
    			cout << "Number of time-steps: " << step << endl;
    		
    		if (step+1 == nstep-1) {  	//At the step before the last write the old final configuration 
    			WriteOldFinal();
    		}
    		
    		if(step%10 == 0){ 				//Measure stuff and write configuration every 10 Monte Carlo steps 
        		Measure();     				//Properties measurement
        		//ConfXYZ(nconf); 			//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        		nconf += 1;   		
	    	}
	    	
	    	ave[iv][iblock] += stima_pot/(double)L;		//Average per block
	     	ave[ik][iblock] += stima_kin/(double)L;	
	     	ave[ie][iblock] += stima_etot/(double)L;	
	     	ave[it][iblock] += stima_temp/(double)L;
	    }
	    	
		for(int k = 0; k < n_props; k++)
			av2[k][iblock] = pow(ave[k][iblock],2);	
  	}
  	
  	for(int i = 0; i < n_props; i++) {
  		for(int j = 0; j < nblocks; j++) {
        	sum_prog = 0;
			sum2_prog = 0;
			for(int l = 0; l < j+1; l++){
				sum_prog += ave[i][l] / (j+1);
				sum2_prog += av2[i][l] / (j+1);
			}
			ave_file[i] << sum_prog << "\t" << Error(sum_prog, sum2_prog, j) << endl;
		}
 	 }
 	for(int i = 0; i < n_props; i++){
		delete[] ave[i];
		delete[] av2[i];
		ave_file[i].close();
  	}
  	
  	delete[] ave;
  	delete[] ave_file;
		 
  	ConfFinal();         			//Write final configuration to restart
  	return 0;
}

void Input(void){ 							//Prepare all stuff for the simulation
	ifstream ReadInput,ReadConf,ReadOld;	
	ofstream scale_factor;
	double ep, ek, pr, et, vir;

	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;

	seed = 1;    		//Set seed for random numbers
	srand(seed); 		//Initialize random number generator
	  
	ReadInput.open("input.dat"); 	//Read input

	ReadInput >> temp;				//Read temperature		
	ReadInput >> npart;				//Read number of particles
	cout << "Number of particles = " << npart << endl;

	ReadInput >> rho;				//Read density
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;		//Compute volume (through density)
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0); 		//Compute edge of the box
	cout << "Edge of the simulation box = " << box << endl;

	ReadInput >> rcut;				//Read cutoff radius of the interaction
	ReadInput >> delta;				//Read delta for time steps
	ReadInput >> nstep;				//Read number of steps
	ReadInput >> iprint;			//Read iprint
	ReadInput >> restart;			//Read wether to restart or not
	
	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl << endl;
	ReadInput.close();

	//Prepare array for measurements, these are indexes
	iv = 0; 						//Potential energy
	ik = 1; 						//Kinetic energy
	ie = 2; 						//Total energy
	it = 3; 						//Temperature
	n_props = 4; 					//Number of observables

	if (restart == 0) { 			//Read initial configuration
		cout << "Start from a new configuration. Read initial configuration from file config.0 " << endl << endl;
	  	ReadConf.open("config.0");
	  	for (int i=0; i < npart; ++i){				//Cycle over each particle	
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = x[i] * box;						//Put every particle inside the box
			y[i] = y[i] * box;
			z[i] = z[i] * box;
	  	}
	  	ReadConf.close();
		
		//Prepare initial velocities
	   	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	   	double sumv[3] = {0.0, 0.0, 0.0};
	   	for (int i=0; i < npart; ++i){					//Choose the velocities randomly	
			vx[i] = rand()/double(RAND_MAX) - 0.5;
		 	vy[i] = rand()/double(RAND_MAX) - 0.5;
		 	vz[i] = rand()/double(RAND_MAX) - 0.5;

		 	sumv[0] += vx[i];							//Sum of the velocities in each direction 
		 	sumv[1] += vy[i];
		 	sumv[2] += vz[i];
	   	}
	   	for (int idim = 0; idim < 3; ++idim) 
	   		sumv[idim] /= (double)npart;
	   
	   	double sumv2 = 0.0, fs;
	  	for (int i=0; i < npart; ++i){					//Subtract the sum of the velocities from the randomly chosen starting velocities in order to have a center of mass velocity equal to zero 
			vx[i] = vx[i] - sumv[0];
		 	vy[i] = vy[i] - sumv[1];
		 	vz[i] = vz[i] - sumv[2];
		 	sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];	//Compute the squared velocities
	   	}
	   	sumv2 /= (double)npart;

	   	scale_factor.open("scale.dat",ios::app);
	   	fs = sqrt(3 * temp/ sumv2);   					// fs = velocity scale factor. Rescales the velocities so that they represent the correct temperatures 
	   	cout << "temp = " << temp << " fs = " << fs << endl;
	   	scale_factor << fs << endl;
	   	scale_factor.close(); 					
	   	
	   	for (int i = 0 ; i < npart; ++i){				//Get the old positions
			vx[i] *= fs;
			vy[i] *= fs;
		 	vz[i] *= fs;

		 	xold[i] = Pbc(x[i] - vx[i] * delta); 		//x(old) = x(actual) - vx(actual)*dt
		 	yold[i] = Pbc(y[i] - vy[i] * delta);
		 	zold[i] = Pbc(z[i] - vz[i] * delta);
	   	}
	} else { 											//Read present and old configuration
		cout << "Restart: Read initial configuration from file old.0 " << endl << endl;
		ReadConf.open("old.0");
		for (int i = 0; i < npart; ++i) {
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = x[i] * box;
			y[i] = y[i] * box;
			z[i] = z[i] * box;
		}
		ReadConf.close();

		cout << "Read old configuration from file old.final " << endl << endl;
		ReadOld.open("old.final");
		for (int i = 0; i < npart; ++i) {
			ReadOld >> xold[i] >> yold[i] >> zold[i];		//Put the particles inside the box
			xold[i] = xold[i] * box;	
			yold[i] = yold[i] * box;
			zold[i] = zold[i] * box;
		}
		ReadOld.close();

		double sumv2 = 0.0, T = 0.0, fs = 0.0;
		for (int i = 0; i < npart; ++i) {
			vx[i] = Pbc(x[i] - xold[i]) / delta;			//Compute the velocities: v = ( x(actual) - x(old) )/dt
			vy[i] = Pbc(y[i] - yold[i]) / delta;
			vz[i] = Pbc(z[i] - zold[i]) / delta;
			sumv2 += vx[i]*vx[i] + vy[i]* vy[i] + vz[i] * vz[i];
		}
		sumv2 /= (double)npart;								//Compute the square velocity per particle
		
		scale_factor.open("scale.dat",ios::app);
	   	fs = sqrt(3 * temp/ sumv2);   					// fs = velocity scale factor. Rescales the velocities so that they represent the correct temperatures 
	   	cout << "temp = " << temp << " fs = " << fs << endl;
	   	scale_factor << fs << endl;
	   	scale_factor.close();
	   	
		for (int i = 0; i < npart; ++i) {
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;
			xold[i] = Pbc(x[i] - vx[i] * delta);
			yold[i] = Pbc(y[i] - vy[i] * delta);
			zold[i] = Pbc(z[i] - vz[i] * delta);
		}
	}
	return;
}

void Move(void){ 												//Move particles with Verlet algorithm
  
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ 									//Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i = 0; i < npart; ++i){ 								//Verlet integration scheme for every particle

  	xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );	//Verlet movement for each direction
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);					//Compute the velocities for each direction
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];												//Save the present coordinates as old
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;												//Save the new coordinates as present
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ 								//Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  							//Distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){											//If the distance between the particles is less than the cutoff radius--> compute the force
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); 	//-Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ 						//Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);	//Open the files in append format
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; 									//Reset observables
  t = 0.0;

  for (int i = 0; i < npart-1; ++i){ 		//Cycle over pairs of particles
    for (int j = i+1; j < npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);							//Compute distance

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       v += vij; 							//Potential energy
     }
    }          
  }

	//Kinetic energy
	for (int i=0; i<npart; ++i) { 
		t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]); 	//Accumulate the kinetic energy for each particle
   }
    stima_pot = v/(double)npart; 								//Potential energy per particle
    stima_kin = t/(double)npart; 								//Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; 				//Temperature: deduce from the kinetic energy
    stima_etot = (t+v)/(double)npart; 							//Total energy per particle
    Epot << stima_pot  << endl;									//Write on files
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    return;
}


void ConfFinal(void){ 												//Write final configuration in "config.final"
  	ofstream WriteConf;

	cout << "Print old configuration to file old.0 " << endl << endl;
	WriteConf.open("old.0");

  	for (int i=0; i<npart; ++i){
    	WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  	}
  	WriteConf.close();
  	return;
}

void WriteOldFinal(void) { //Write final configuration
	ofstream WriteOld;
	cout << "Print final configuration to file old.final " << endl << endl;
  	WriteOld.open("old.final");
	for (int i = 0; i < npart; ++i)
		WriteOld << x[i] / box << " " << y[i] / box << " " << z[i] / box << endl;
	WriteOld.close();
return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  	ofstream WriteXYZ;
  	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  	WriteXYZ << npart << endl;
  	WriteXYZ << "This is only a comment!" << endl;
  	for (int i=0; i<npart; ++i){
    	WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  	}
  	WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r/box);
}

double Error (double ave, double ave2, int i) {
	if (i == 0)
		return 0;
	else 
		return sqrt((ave2-pow(ave,2))/i);
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

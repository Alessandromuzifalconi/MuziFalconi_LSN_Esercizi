/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie, igofr;
double stima_pot, stima_kin, stima_etot, stima_temp, err_gdir;
double bin_size,nbins;
double restart;
double walker[m_props];
// averages
double acc,att;
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
const int nblocks = 100;
double err_pot;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];		//Position for every particle both at present time step and previous time step (old)
double vx[m_part],vy[m_part],vz[m_part];											//Speed for every particle

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;		//Energy and temperature self explicatory, vol = volume, rho = density, box = side of the box containing the particles, rcut = cutoff radius

// simulation
int nstep, nblk, iprint, seed;		//nstep = number of steps fo the simulation (10000), iprint = print the step every iprint (1000) number of steps, seed for generating the random numbers
double delta;	

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void WriteOldFinal(void); 							//Write old final configuration
void Reset(int);
void Accumulate(void);
void Averages(int);
double Error(double, double, int);					//Error for the blocking method
 
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

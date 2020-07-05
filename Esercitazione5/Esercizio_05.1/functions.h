Random rnd; 
int seed[4];
int p1, p2;
const int M = 1000000; 					//Number of steps
const int N = 200; 						//Number of blocks
int L = M/N; 							//Number of elements in each block
const int nsave = 2;					//One every nsave configurations is saved for showing the config. 
int n_acc = 0, n_rej = 0; 				//Number of accepted and rejected steps
Position x1; 							//Starting position
Position x2; 							//Second posiion
double r_max; 							//Maximum size of the step 
double alpha;							//Variable for implementing the minimum
double r[M];

double *ave = new double[N] {0}; 		//Variables for the blocking method
double *ave2 = new double[N] {0};			
double *sum_prog = new double [N] {0};
double *sum2_prog = new double [N] {0};
double *err_prog = new double [N] {0};

double Psi100 (Position); 				//Returns the probability density of point r for the wavefunction Psi 100. Returns the square modulus of Psi100
double Psi210 (Position); 				//Returns the probability density of point r for the wavefunction Psi 210 Returns the square modulus of Psi210

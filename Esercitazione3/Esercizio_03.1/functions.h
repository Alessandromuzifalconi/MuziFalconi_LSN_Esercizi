Random rnd;
int seed[4];
int p1, p2;
const int N = 50000;					//Number of times we repeat the simulation
const int steps = 100;					//Number of steps for each simulation = number of blocks
int L = N/steps;						//Values for each block

const double T = 1.; 					//Delivery time
const double S0 = 100.;					//Price at t = 0
const double K = 100.;					//Strike price
const double r = 0.1;					//Interest rate
const double sigma = 0.25;				//Volatility
double ST = 0;							//Value of S(T)
double z;

double *S = new double [steps];			//Price at each time step
double *t = new double [steps]; 		//Time steps

double *C = new double [L] {0};				//Call option price, one for each step in the block
double *P = new double [L] {0};				//Put option price, one for each step in the block

double *sum_prog = new double [steps] {0};	//Variables for the blocking method
double *sum2_prog = new double [steps] {0};
double *err_prog = new double [steps] {0};
double *ave = new double [steps] {0};
double *ave2 = new double [steps] {0};

double *sump_prog = new double [steps] {0};	//Variables for the blocking method: put asset
double *sum2p_prog = new double [steps] {0};
double *errp_prog = new double [steps] {0};
double *avep = new double [steps] {0};
double *ave2p = new double [steps] {0};

void Random_Initializer();					//Initialize pseudo random number generator

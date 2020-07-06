Random rnd;
int seed[4];
int p1;
int p2;

const int N = 32; 					//Number of cities
double cities[N][2];				//Positions of the cities
int method = 1; 					//If method is 0 we generate the cities on a circumference, otherwise we generate them inside a square.

const int m = 100;					//Number of temperature visited
const int n = 10000;				//Steps for each temperature
double beta = 0;					//Inverse temperature
double Tf = 0.01;					//Final temperature

int n_acc = 0, n_rej = 0; 			//Number of accepted and rejected steps
double alpha;						//Variable for implementing the minimum
int* seq_old = new int[N] {0}; 		//Sequence for initializing
int i_s = 1, j_s = 1;				//Indexes


void RandomInitializer();
void Cities_generator (int, int); 			//Generates the position of each city
int Check (int*); 							//Checks if the individual satisfies the right boundaries: returns 1 if the individual is not ok, returns 0 otherwise
int PBC(int); 								//Periodic boundary condition
double L (int*); 							//Computes the distance travelled by the salesman for a given sequence of cities
void Exchange (int*); 						//Permutes two cities in an individual

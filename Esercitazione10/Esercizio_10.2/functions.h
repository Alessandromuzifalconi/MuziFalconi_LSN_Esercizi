Random rnd;
int seed[4];
int p1;
int p2;

const int N = 32; 				//Number of cities
int P = 500; 					//Size of the population
double cities[N][2];			//Position of the cities
int method = 1; 				//If method is 0 we generate the cities on a circumference, otherwise we generate them inside a square.
const int Generations = 1000; 	//Number of gen	erations
const int Nmigr = 100;			//Number of migrations
int Mig = Generations/Nmigr;	//After ho	w many generations there is a migration

int sent[N] {0};				//Sequence to be sent
int received[N] {0};			//Sequence to be received
int n = 1;						//Just counts the generation

double* lenghts = new double [P] {0};		//Distance travelled by each individual
int* order = new int [P];					//Order of the individuals
int r = 0;									//Parameter for determining the migration

//Mutation probabilities
double P_PairExchange = 0.1;
double P_Shift = 0.1;
double P_Permutation = 0.1;
double P_Inversion = 0.1;
double P_Crossing = 0.5;

void RandomInitializer(int);
void Cities_generator (int method, int N); 			//Generates the position of each city
void Population_generator (Individual*); 			//Generates the starting population (each individual has a random sequence)
int Check (Individual*, int); 						//Checks if the individual satisfies the right boundaries: returns 1 if the individual is not ok, returns 0 otherwise
int PBC(int); 										//Periodic boundary condition
double L (int* sequence); 							//Computes the distance travelled by the salesman for a given sequence of cities
void L_Calculator (Individual*, double*);			//Computes the distance for every individual
void Quicksort(double*, int*, int, int);			//Orders from minimum to maximum
void Order_Individuals (Individual*, int*);			//Orders individuals from best to worst
int Selection (int);								//Select an individual for crossing over
int Selection_Mutation (int);						//Select an individual for single individual-mutations
void Mutation_PairExchange (Individual*, int); 		//Permutes two cities in an individual
void Mutation_Shift(Individual*, int);				//Shift cities in an individual
void Mutation_Permutation(Individual*, int);		//Permutes more than two cities in an individual
void Mutation_Inversion(Individual*, int);			//Invert order of cities in an individual
void CrossingOver(Individual*, int, int); 			//Do crossing over between two individuals
void Mutate (Individual*);							//Function for deciding whether to mutate or not, and what mutation to make

//Migrations
double Determine_Mig ();							//Determine who migrates where
void Migration_A (int*, int*, int);					//0 <---> 1, 2 <---> 3 
void Migration_B (int*, int*, int);					//0 <---> 2, 1 <---> 3 
void Migration_C (int*, int*, int);					//0 <---> 3, 1 <---> 2 

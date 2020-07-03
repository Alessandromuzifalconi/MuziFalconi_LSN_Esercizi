Random rnd;
int seed[4];
int p1, p2;
int M = 10000;										//Number of experiments 
double *S1 = new double [M]; 						//Array to be filled with S1 values
double *S2 = new double [M];						//Array to be filled with S2 values
double *S10 = new double [M];						//Array to be filled with S10 values
double *S100 = new double [M];						//Array to be filled with S100 values
double x;											//Variable where I save the random numbers when needed
int N;												//Number of sums (N=1,2,10,100)


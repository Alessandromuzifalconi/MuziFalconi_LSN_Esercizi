Random rnd;
int seed[4];
int p1, p2;
Position pos(0,0,0); 					//The random walk starts at the origin (0,0,0)
int N = 100000; 						//Number of times we repeat the random walk
int steps = 100; 						//Number of steps for each random walk
double a = 1; 							//Size of the step
double x,theta,phi;						//x is just an helpful variable for saving things when needed, theta and phi are angles for the continuos random walk
double *r = new double[steps];			//Distance from the origin 

double *sum_prog = new double [steps] {0};	//Variables for the blocking method
double *sum2_prog = new double [steps] {0};
double *err_prog = new double [steps] {0};
double *sums = new double [steps] {0};
double *sums2 = new double [steps]  {0};

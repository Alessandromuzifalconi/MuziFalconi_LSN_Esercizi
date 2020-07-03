Random rnd;
int seed[4];
int p1, p2;
int M = 100000; 	//umber of throws
int N = 200; 		//Number of blocks for the blocking method
int K = M/N; 		//Number of throws in each block
double x;			//x-coordinate of the left tip of the needle (x-axis is orthogonal to the lines that divide the area)
double l;			//x projection of the needle
double theta;		//Angle between the needle and the x-axis
double d = 1.;		//Spacing between the lines 
double L = 0.7;		//Length of the needle
int hits = 0;		//Number of hits 

double *pi = new double[N] {0}; 		//Variables for the blocking method 
double *pi2 = new double[N] {0};			
double *sum_prog = new double [N] {0};
double *sum2_prog = new double [N] {0};
double *err_prog = new double [N] {0};
	

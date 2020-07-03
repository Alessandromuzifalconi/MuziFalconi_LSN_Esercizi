Random rnd;
int seed[4];
int p1, p2;
int M = 40000; 		//Number of throws (i.e. pseudo numbers generated)
int N = 100; 		//Number of blocks for the blocking method
int L = M/N; 		//Number of throws in each block
int m = 100; 		//Number of sub-intervals in which we divide [0,1] for the Xi2 test
int n = 10000; 		//Number of attempted throws in each sub-interval for implementig the Xi2 test
int counter; 		//Number of throws that actually fall in the sub-interval
double x;			//Just a variable where I save the random number (when I need to)


double *ave = new double [N] {0};			//Variables for the blocking method for the mean value
double *ave2 = new double [N] {0};
double *sum_prog = new double [N] {0};
double *sum2_prog = new double [N] {0};
double *err_prog = new double [N] {0};
double sum = 0, sum_var = 0;
	
double *ave_var = new double [N] {0}; 		//Variables for the blocking method for the variance 
double *ave2_var = new double [N] {0};
double *sum_prog_var = new double [N] {0};
double *sum2_prog_var = new double [N] {0};
double *err_prog_var = new double [N] {0};
	
double *count = new double [m] {0}; 			//Number of throws in each sub interval 
double *Xi2_block = new double [N] {0}; 		//Xi2 for each block
double *Xi2_block2 = new double [N] {0}; 		//squared Xi2 for each block
	

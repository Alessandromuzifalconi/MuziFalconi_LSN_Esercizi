	//Pseudo-random number generator
	Random rnd;
	int seed[4];
	int p1, p2;
	
	//Averages
	const int M = 200000; 							//Number of steps
	const int N = 500; 							//Number of blocks
	int L = M/N; 									//Number of elements in each block
	
	//Metropolis
	double n_acc = 0.0, n_rej = 0.0; 				//Number of accepted and rejected steps
	double x_old, x_new;
	const double x_max = 2.75; 						//Maximum size of the step
	double alpha;									//Variable for evaluating the minimum
	
	double x [M];									//Array with all the positions distributed as |Psi|^2
	double mu, sigma;								//Parameters for optimization
	
	//Variables for the blocking method: uniform
	double *ave = new double[N] {0}; 		
	double *ave2 = new double[N] {0};			
	double *sum_prog = new double [N] {0};
	double *sum2_prog = new double [N] {0};
	double *err_prog = new double [N] {0};
	double sum = 0;
	int k;
	
	//Functions
	double Error (double* ave, double* ave2, int i); 	//Returns statistical uncertainty	
	double Psi (double x, double sigma, double mu); 	//Returns the value of the approximated wavefunction evaluated in x (1D)
	double V (double x);								//Returns the potential evaluated in x
	double HPsi (double x, double sigma, double mu); 	//Returns H applied to Psi

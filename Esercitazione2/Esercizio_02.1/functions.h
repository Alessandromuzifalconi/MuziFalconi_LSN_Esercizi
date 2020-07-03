Random rnd;
int seed[4];
int p1, p2;
int M = 100000; 		//Number of throws
int N = 100; 			//Number of blocks for the blocking method
int L = M/N; 			//Number of throws in each block
double x;				//Just a variable where I save the random number (when I need to)

double *ave = new double[N] {0}; 		//Variables for the blocking method: uniform sampling
double *ave2 = new double[N] {0};			
double *sum_prog = new double [N] {0};
double *sum2_prog = new double [N] {0};
double *err_prog = new double [N] {0};

double *ave_imp = new double[N] {0}; 		//Variables for the blocking method: importance sampling
double *ave2_imp = new double[N] {0};			
double *sum_prog_imp = new double [N] {0};
double *sum2_prog_imp = new double [N] {0};
double *err_prog_imp = new double [N] {0};
double sum = 0, sum_imp = 0;

double importance_sampl (); 			//Generates a random variable distributed as 2(1-x) for importance sampling

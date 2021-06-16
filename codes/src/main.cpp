#include <omp.h>

#include "Para.h"
#include "Statis.h"
#include "Solver.h"

using namespace std;


void Config(int n);


/***** computation mode *****/
# define AUXIMAIN // DEFAULT //

#ifdef DEFAULT

int main()
{
	Solver slv("");

	Config(slv.para.nthrds);

	while (slv.step < slv.para.Nt) {

		slv.Evolve();
		slv.Output();

		if (slv.step % slv.para.nprint == 0) {
			// check modifications to XINDAT every nprint steps
			slv.para.checkPara(slv.step);
			Config(slv.para.nthrds);
		}
	}

	cout << "\nComputation finished!" << endl;
}

#endif
#ifdef AUXIMAIN

int main()
{
	Solver slv0("base/");
	Solver slv1("test/");

	Config(slv1.para.nthrds);

	while (slv1.step < slv1.para.Nt) {

		while (slv0.time - slv1.time < slv1.para.dt - INFTSM) {
			slv0.Evolve();
			slv0.Output();
		}

		slv1.Evolve(slv0);
		slv1.Output();

		if (slv1.step % slv1.para.nprint == 0) {
			slv0.para.checkPara(slv0.step);
			slv1.para.checkPara(slv1.step);
			Config(slv1.para.nthrds);
		}
	}

	cout << "\nComputation finished!" << endl;
}

#endif

// notes of parameters used in constant flowrate computations
// slv0.mpg[0] = -2.900e-3; slv1.mpg[0] = -2.9000e-3; // Re 10000, Ret 540
// slv0.mpg[0] = -1.736e-3; slv1.mpg[0] = -2.5000e-3; // Re 20000, Ret 1000
// slv0.mpg[0] = -0.865e-3; slv1.mpg[0] = -2.1140e-3; // Re 43500, Ret 2000


void Config(int n)
{
	static int nthrds = 0;

	int tempn = nthrds;
	int nprocs = omp_get_num_procs();
	// int N = para.Nx * para.Ny * para.Nz;

	if (n > 0)
		nthrds = (n > nprocs ? nprocs : n);
	// else {
	// 	// automatically decide number of OpenMP threads based on grid number
	// 	nthrds = log2(N) - 16;
	// 	nthrds = ( nthrds < 1 ? 1 : (
	// 			   nthrds > 32? 32: (
	// 			   nthrds > nprocs ? nprocs : nthrds )));
	// }
	if (tempn != nthrds) {
		omp_set_num_threads(nthrds);
		cout << endl << "number of OpenMP threads: " << nthrds << endl;
	}
}


/***** time test example *****/
// #include <sys/time.h>

// struct timeval *time0 = new struct timeval;
// struct timeval *time1 = new struct timeval;
// long duration = 0;

// gettimeofday(time0, NULL);
// ///// codes to be timed
// gettimeofday(time1, NULL);

// duration = 1e6 * (time1->tv_sec - time0->tv_sec) + (time1->tv_usec - time0->tv_usec);
/*****************************/







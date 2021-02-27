#include <omp.h>

#include "Para.h"
#include "Statis.h"
#include "Solver.h"

using namespace std;


/***** computation mode *****/
# define DEFAULT // AUXIMAIN //


void Config(int n, const Para &para);
int Initiate(Solver &solver, const Para &para);
void Output(Para &para, Solver &solver, int tstep);


#ifdef DEFAULT

int main()
{
	Para para("");

	Geometry_prdz geo(para.Nx,para.Ny,para.Nz, para.Lx,para.Ly,para.Lz);

	geo.InitMesh(para.dy_min);
	geo.InitInterval();
	geo.InitWaveNumber();
	geo.InitIndices();

	para.showPara();
	geo.WriteMesh(para.statpath);
	geo.WriteMeshY(para.statpath);
	
	const Mesh mesh(geo);
	Solver solver(mesh);

	// computation begins
	int tstep = Initiate(solver, para);

	solver.set_mpg(0, 0, 0);
	
	if (tstep == 0) Output(para, solver, tstep);

	// main loop
	while (tstep++ < para.Nt) {
		solver.evolve(para.Re, para.dt, para.bftype);
		Output(para, solver, tstep);
	}

	cout << "\nComputation finished!" << endl;
}

#endif
#ifdef AUXIMAIN

int main()
{
	Para para0("base/");
	Para para1("test/");

	Geometry_prdxz geo0(para0.Nx,para0.Ny,para0.Nz, para0.Lx,para0.Ly,para0.Lz);
	Geometry_prdxz geo1(para1.Nx,para1.Ny,para1.Nz, para1.Lx,para1.Ly,para1.Lz);

	geo0.InitMesh(para0.dy_min);
	geo1.InitMesh(para1.dy_min);

	// geo1.AlignBoundaryYc(geo0);

	geo0.InitInterval();   geo1.InitInterval();
	geo0.InitWaveNumber(); geo1.InitWaveNumber();
	geo0.InitIndices();    geo1.InitIndices();

	para0.showPara(); geo0.WriteMeshY(para0.statpath); geo0.WriteMesh(para0.statpath);
	para1.showPara(); geo1.WriteMeshY(para1.statpath); geo1.WriteMesh(para1.statpath);

	const Mesh mesh0(geo0);
	const Mesh mesh1(geo1);
	Solver solver0(mesh0);
	Solver solver1(mesh1);

	// computation begins
	int tstep0 = Initiate(solver0, para0);
	int tstep1 = Initiate(solver1, para1);


	solver0.set_mpg(-2.9e-3,   0, 0); solver1.set_mpg(-2.9e-3,    0, 0); // Ret 540
	solver0.set_mpg(-1.736e-3, 0, 0); solver1.set_mpg(-2.5e-3,    0, 0); // Ret 1000
	solver0.set_mpg(-8.65e-4,  0, 0); solver1.set_mpg(-2.114e-3,  0, 0); // Ret 2000
	solver0.set_mpg(-8.65e-4,  0, 0); solver1.set_mpg(-1.9753e-3, 0, 0); // Ret 4000 (MFU 2000)
	solver0.set_mpg(-1.,       0, 0); solver1.set_mpg(-1.,        0, 0); // Ret specified by XINDAT

	
	if (tstep0 == 0) Output(para0, solver0, tstep0);
	if (tstep1 == 0) Output(para1, solver1, tstep1);

	// main loop
	while (tstep1++ < para1.Nt) {

		while (solver0.get_time() - solver1.get_time() < para1.dt - INFTSM) {
			tstep0 ++;
			solver0.evolve(para0.Re, para0.dt, para0.bftype);
			Output(para0, solver0, tstep0);
		}

		solver1.evolve(para1.Re, para1.dt, para1.bftype, solver0, para0.Re);
		Output(para1, solver1, tstep1);
	}

	cout << "\nComputation finished!" << endl;
}

#endif


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


int Initiate(Solver &solver, const Para &para)
{
	int tstep = 0;
	Config(para.nthrds, para);
	solver.set_time(0);

	if (para.nread == 0) {
		solver.init(para.inener);
		cout << endl << "Flow fields initiated from laminar." << endl;
	}
	else if (! strcmp(Para(para.inpath).statpath, para.statpath)) {
		tstep = para.nread;
		solver.set_time(Statis::GetLogTime(para.statpath, tstep, solver.get_mpg()));
		solver.get_fld().ReadField(para.fieldpath, tstep, "");
		cout << endl << "Continue from step " << tstep << ", time " << solver.get_time() << endl;
	}
	else {
		Para para0(para.inpath);
		Geometry geo0(para0.Nx, para0.Ny, para0.Nz, para0.Lx, para0.Ly, para0.Lz);
		geo0.InitMesh(0, para0.statpath);
		// geo0.InitIndices();

		solver.init(Flow(Mesh(geo0)).ReadField(para0.fieldpath, para.nread, ""));

		cout << endl << "Flow fields initiated from existing fields with parameters:" << endl;
		para0.showPara();
	}

	return tstep;
}

void Output(Para &para, Solver &solver, int tstep)
{	
	if (tstep % para.nwrite == 0) {
		solver.get_fld ().WriteField(para.fieldpath, tstep, "");
		solver.get_fldh().WriteField(para.fieldpath, tstep, "T");
		solver.get_fld().WriteTecplot(para.fieldpath, 0, solver.get_time());
		cout << "Files successfully written for step " << tstep << endl;
	}
	if (tstep % para.nprint == 0) {

		Statis statis(solver.ms);

		statis.Check(solver.get_fld(), solver.get_vis(), para.Re, para.dt);
		statis.WriteProfile(para.statpath);
		statis.WriteLogfile(para.statpath, tstep, solver.get_time(), solver.get_mpg());

		para.checkPara(tstep);
		
		Config(para.nthrds, para);
	}
}

void Config(int n, const Para &para)
{
	static int nthrds = 0;

	int tempn = nthrds;
	int nprocs = omp_get_num_procs();
	int N = para.Nx * para.Ny * para.Nz;

	if (n > 0)
		nthrds = (n > nprocs ? nprocs : n);
	else {
		// automatically decide number of OpenMP threads based on grid number
		nthrds = log2(N) - 16;
		nthrds = ( nthrds < 1 ? 1 : (
				   nthrds > 32? 32: (
				   nthrds > nprocs ? nprocs : nthrds )));
	}
	if (tempn != nthrds) {
		omp_set_num_threads(nthrds);
		cout << endl << "number of OpenMP threads: " << nthrds << endl;
	}
}







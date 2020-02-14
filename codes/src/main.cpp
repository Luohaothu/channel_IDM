# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <ctime>

# include "Para.h"
# include "Statis.h"
# include "Solver.h"

using namespace std;




void initiate(Solver &solver, int &tstep, double &time, const Para &para)
{
	tstep = 0;
	time = 0;
	solver.config(para.nthrds);

	if (para.nread == 0) {
		solver.FLD.initrand(para.inener);
		cout << endl << "Flow fields initiated from laminar." << endl;
	}
	else if (! strcmp(Para(para.inpath).statpath, para.statpath)) {
		tstep = para.nread;
		time = Statis::getLogtime(para.statpath, tstep);
		solver.FLD.readField(para.fieldpath, tstep, "");
		cout << endl << "Continue from step " << tstep << ", time " << time << endl;
	}
	else {
		Para para0(para.inpath);
		Mesh mesh0(para0.Nx, para0.Ny, para0.Nz, para0.Lx, para0.Ly, para0.Lz, 0, para0.statpath);
		solver.FLD.initfrom(Feld(mesh0).readField(para0.fieldpath, para.nread, ""));
		cout << endl << "Flow fields initiated from existing fields with parameters:" << endl;
		para0.showPara();
		mesh0.freeall();
	}
}

void output(Statis &statis, Para &para, Solver &solver, int tstep, double time)
{
	if (tstep % para.nwrite == 0) {
		solver.FLD.writeField(para.fieldpath, tstep, "");
		solver.FLDH.writeField(para.fieldpath, tstep, "T");
		// solver.FLD.writeTecplot(para.fieldpath, tstep, time);
		cout << "Files successfully written for step " << tstep << endl;
	}
	if (tstep % para.nprint == 0) {
		statis.check(solver.FLD, solver.VIS.S, para.Re, para.dt);
		statis.writeProfile(para.statpath);
		statis.writeLogfile(para.statpath, tstep, time, solver.mpg);

		para.checkPara(tstep);
		solver.config(para.nthrds);
	}
}


# define DATAASM

# ifndef DATAASM

int main()
{
	Para para("");

	Mesh mesh(para.Nx,para.Ny,para.Nz, para.Lx,para.Ly,para.Lz, para.dy_min);
	Mesh bmesh(mesh.Nx,2-1,mesh.Nz, mesh.Lx,mesh.Ly,mesh.Lz);

	bmesh.y[0] = 0; bmesh.yc[0] = mesh.y[1];
	bmesh.y[1] = 0; bmesh.yc[1] = mesh.y[mesh.Ny];

	Solver solver(mesh, bmesh); Statis statis(mesh);

	int tstep = 0; double time = 0;

	// computation begins
	para.showPara(); mesh.writeYmesh(para.statpath);

	initiate(solver, tstep, time, para);
	if (tstep == 0) output(statis, para, solver, tstep, time);

	// main loop
	while (tstep++ < para.Nt) {
		time += para.dt;
		solver.evolve(para.Re, para.dt, para.bftype);
		output(statis, para, solver, tstep, time);
	}

	cout << "\nComputation finished!" << endl;
}

# else

int main()
{
	Para para0("experiment/"), para1("dataasm/");

	Mesh mesh0(para0.Nx,para0.Ny,para0.Nz, para0.Lx,para0.Ly,para0.Lz, para0.dy_min);
	Mesh mesh1(para1.Nx,para1.Ny,para1.Nz, para1.Lx,para1.Ly,para1.Lz, para1.dy_min);

	Mesh bmesh0(mesh0.Nx,2-1,mesh0.Nz, mesh0.Lx,mesh0.Ly,mesh0.Lz);
	Mesh bmesh1(mesh1.Nx,2-1,mesh1.Nz, mesh1.Lx,mesh1.Ly,mesh1.Lz);

	bmesh0.y[0] = 0; bmesh0.yc[0] = mesh0.y[1];
	bmesh0.y[1] = 0; bmesh0.yc[1] = mesh0.y[mesh0.Ny];
	bmesh1.y[0] = 0; bmesh1.yc[0] = mesh1.y[1];
	bmesh1.y[1] = 0; bmesh1.yc[1] = mesh1.y[mesh1.Ny];

	Solver solver0(mesh0, bmesh0); Statis statis0(mesh0);
	Solver solver1(mesh1, bmesh1); Statis statis1(mesh1);

	int tstep0 = 0; double time0 = 0;
	int tstep1 = 0; double time1 = 0;

	DA da(mesh1);

	// computation begins
	para0.showPara(); mesh0.writeYmesh(para0.statpath);
	para1.showPara(); mesh1.writeYmesh(para1.statpath);

	initiate(solver0, tstep0, time0, para0);
	initiate(solver1, tstep1, time1, para1);
	if (tstep0 == 0) output(statis0, para0, solver0, tstep0, time0);
	if (tstep1 == 0) output(statis1, para1, solver1, tstep1, time1);

	// main loop
	while (tstep1++ < para1.Nt) {
		time1 += para1.dt;

		while (time1 - time0 > INFTSM) {
			tstep0 ++;
			time0 += para0.dt;
			solver0.evolve(para0.Re, para0.dt, para0.bftype);
			output(statis0, para0, solver0, tstep0, time0);
		}

		/* FLOWS COMMUNICATION: introduce experiment data */
		if (da.getExp(time1, solver0.FLD.V))
			solver1.evolve_da(para1.Re, para1.dt, para1.bftype, da);
		else
			solver1.evolve   (para1.Re, para1.dt, para1.bftype);

		output(statis1, para1, solver1, tstep1, time1);
	}

	cout << "\nComputation finished!" << endl;
}

# endif




/***** time test example *****/
// # include <sys/time.h>

// struct timeval *time0 = new struct timeval;
// struct timeval *time1 = new struct timeval;
// long duration = 0;

// gettimeofday(time0, NULL);
// ///// codes to be timed
// gettimeofday(time1, NULL);

// duration = 1e6 * (time1->tv_sec - time0->tv_sec) + (time1->tv_usec - time0->tv_usec);
/*****************************/





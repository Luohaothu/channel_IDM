# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <ctime>
# include <sys/stat.h>

# include "Para.h"
# include "Statis.h"
# include "Solver.h"

using namespace std;


class Flow
{
	public:
		Para para;
		Mesh mesh;
		Field field;
		Statis statis;
		Solver solver;
		int tstep;
		double time;

		Flow(char path[1024]):
		para(path),
		mesh(para.Nx,para.Ny,para.Nz, para.Lx,para.Ly,para.Lz, para.dy_min),
		field(mesh),
		statis(mesh),
		solver(mesh, field)
		{
			tstep = 0;
			time = 0;
			XINDAT_modify_time = 0;
			strcpy(workpath, path);

			para.showPara();
			mesh.writeYmesh(para.statpath);

			this->initiate();
			this->input(); // for solver config
			if (tstep == 0)	this->output();
		};

		Field& solve(double target_time)
		{
			while (target_time - time > 1e-10) {
				tstep ++;
				time += para.dt;
				solver.getup(para.Re, para.dt, para.bftype);
				this->output();
			}
			return field;
		};

		void output()
		{
			if (tstep % para.nwrite == 0) {
				field.writeField(para.fieldpath, tstep);
				cout << "Files successfully written for step " << tstep << endl;
			}
			if (tstep % para.nprint == 0) {
				statis.check(field.U, field.P, field.NU, para.Re, para.dt);
				statis.writeProfile(para.statpath);
				statis.writeLogfile(para.statpath, tstep, time, field.mpg);

				this->input();
			}
		};

	private:
		time_t XINDAT_modify_time;
		char workpath[1024];
		
		void initiate()
		{
			if (para.nread == 0) {
				field.initField(para.inener);
				cout << endl << "Flow fields initiated from laminar." << endl;
			}
			else if (! strcmp(Para(para.inpath).statpath, para.statpath)) {
				tstep = para.nread;
				time = statis.getLogtime(para.statpath, tstep);
				field.readField(para.fieldpath, tstep);
				cout << endl << "Continue from step " << tstep << ", time " << time << endl;
			}
			else {
				Para para0(para.inpath);
				Mesh mesh0(para0.Nx, para0.Ny, para0.Nz, para0.Lx, para0.Ly, para0.Lz, 0, para0.statpath);

				field.initField(Field(mesh0).readField(para0.fieldpath, para.nread));

				cout << endl << "Flow fields initiated from existing fields with parameters:" << endl;
				para0.showPara();
				mesh0.freeall();
			}
		};

		void input()
		{
			char filename[1024]; strcat(strcpy(filename, workpath), "XINDAT");
			struct stat buf; stat(filename, &buf);
			time_t tempt = XINDAT_modify_time;

			if (tempt < (XINDAT_modify_time = buf.st_mtime)) {

				para.readPara(workpath);

				if ((int)tempt != 0) {
					cout << endl << filename << " updated at step " << tstep << " as:" << endl;
					para.showPara();
				}
			}

			solver.config(para.nthrds);
		};
};


// # define CONTROL_TYP 4

# ifndef CONTROL_TYP
/***** wall-boundary computation *****/

int main()
{
	Flow flow("");
	while (flow.tstep < flow.para.Nt)
		flow.solve(flow.time + flow.para.dt);
	cout << "Computation finished!" << endl;
}

# else
/***** off-wall boundary computation *****/

int main()
{
	Flow flow("offwall/"), flow0("bcgenerator/");

	while (flow.tstep++ < flow.para.Nt) {
		flow.time += flow.para.dt;

		flow.field.bcond(flow0.solve(flow.time).U);			// set boundary conditions
		flow.field.getnu(flow.para.Re, flow.para.bftype);	// get the viscosity field
		flow.field.getup(flow.para.dt, flow.para.nthrds);	// time evolution
		if (flow.para.bftype == 1) flow.field.removeSpanMean();	// for MFU
		
		/* flow rate modification for off-wall calculation */
		if (2.0 - flow.para.Ly > 1e-10) {

// type 4: keep the pressure gradients equal to that of the full-sized channel
			double dmpg1 = flow0.field.mpg[0] - flow.field.mpg[0];
			double dmpg3 = flow0.field.mpg[2] - flow.field.mpg[2];
			flow.field.mpg[0] += dmpg1;
			flow.field.mpg[2] += dmpg3;
			for (int j=1; j<flow.para.Ny; j++) {
				flow.field.U.com1.layerAdd(- flow.para.dt * dmpg1, j);
				flow.field.U.com3.layerAdd(- flow.para.dt * dmpg3, j);
			}

		}


		flow.output();
		flow.input ();
	}
}

# endif
/***** ***** *****/





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





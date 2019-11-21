# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <ctime>
# include <sys/stat.h>

# include "Para.h"
# include "Field.h"
# include "Statis.h"

using namespace std;


class FluidSolver
{
	public:
		Para para;
		Mesh mesh;
		Field field;
		Statis stas;
		int tstep;
		double time;

		FluidSolver(char path[1024]):
			para(path),
			mesh(para.Nx,para.Ny,para.Nz, para.Lx,para.Ly,para.Lz, para.dy_min),
			field(mesh),
			stas(mesh),
			tstep(0),
			time(0)
			{
				para.showPara();
				mesh.writeYmesh(para.statpath);
				this->initiate();
				if (tstep == 0)	this->output();
				
				XINDAT_modify_time = 0;
				strcpy(workpath, path);
			};

		Field& solve(double target_time) {
			while (target_time - time > 1e-10) {
				tstep ++;
				time += para.dt;

				field.bcond(tstep);					// set boundary conditions
				field.getnu(para.Re, para.bftype);	// get the viscosity field
				field.getup(para.dt, para.nthrds);	// time evolution
				if (para.bftype == 1) field.removeSpanMean();	// for MFU
				
				this->output();
				this->input ();
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
				stas.check(field.U, field.P, field.NU, field.mpg, para.Re, para.dt);
				stas.writeProfile(para.statpath);
				stas.writeLogfile(para.statpath, tstep, time);
			}
		};

		void input()
		{
			if (tstep % para.nprint == 0) {
				char filename[1024];
				struct stat buf;
				stat(strcat(strcpy(filename, workpath), "XINDAT"), &buf);

				if (buf.st_mtime > XINDAT_modify_time && (int)XINDAT_modify_time != 0) {
					para.readPara(workpath);
					cout << endl << filename << " updated at step " << tstep << " as:" << endl;
					para.showPara();
				}
				XINDAT_modify_time = buf.st_mtime;
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
				time = stas.getLogtime(para.statpath, tstep);
				field.readField(para.fieldpath, tstep);
				cout << endl << "Continue from step " << tstep << ", time " << time << endl;
			}
			else {
				Para para0(para.inpath);
				Mesh mesh0(para0.Nx, para0.Ny, para0.Nz, para0.Lx, para0.Ly, para0.Lz, 0, para0.statpath);
				Field field0(mesh0);
				field0.readField(para0.fieldpath, para.nread);
				field.initField(field0);
				mesh0.freeall();
				cout << endl << "Flow fields initiated from existing fields with parameters:" << endl;
				para0.showPara();
			}
		};
};

int main()
{
	FluidSolver flow("");
	while (flow.tstep++ < flow.para.Nt) {
		flow.time += flow.para.dt;

		flow.field.bcond(flow.tstep);						// set boundary conditions
		flow.field.getnu(flow.para.Re, flow.para.bftype);	// get the viscosity field
		flow.field.getup(flow.para.dt, flow.para.nthrds);	// time evolution
		if (flow.para.bftype == 1) flow.field.removeSpanMean();	// for MFU
		
		flow.output();
		flow.input ();
	}
}






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





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


# define CONTROL_TYP 2

int main()
{
	FluidSolver flow("offwall/"), flow0("bcgenerator/");

	while (flow.tstep++ < flow.para.Nt) {
		flow.time += flow.para.dt;

		flow.field.bcond(flow0.solve(flow.time).U);			// set boundary conditions
		flow.field.getnu(flow.para.Re, flow.para.bftype);	// get the viscosity field
		flow.field.getup(flow.para.dt, flow.para.nthrds);	// time evolution
		if (flow.para.bftype == 1) flow.field.removeSpanMean();	// for MFU
		
		/* flow rate modification for off-wall calculation */
		if (2.0 - flow.para.Ly > 1e-10) {

# if CONTROL_TYP == 1
// type 1: fix the streamwise flowrate to be constant and spanwise 0
			double dmpg = (flow.field.U.com1.bulkMeanU() - 1.01876235787727) / flow.para.dt;
			flow.field.mpg[0] += dmpg;
			for (int j=1; j<flow.para.Ny; j++) flow.field.U.com1.layerAdd(- flow.para.dt * dmpg, j);

# elif CONTROL_TYP == 2
// type 2: dynamically adjust mean pressure gradients to keep the flowrates identical to that of the full-sized channel
			int jb1 = (flow0.para.Ny - flow.para.Ny) / 2, jb2 = flow0.para.Ny - jb1;
			double um = 0, wm = 0;

			for (int j=jb1+1; j<jb2; j++) {
				um += flow0.field.U.com1.layerMean(j) * flow0.mesh.dy[j] / flow.mesh.Ly;
				wm += flow0.field.U.com3.layerMean(j) * flow0.mesh.dy[j] / flow.mesh.Ly;
			}
			double dmpg1 = (flow.field.U.com1.bulkMeanU() - um) / flow.para.dt;
			double dmpg3 = (flow.field.U.com3.bulkMeanU() - wm) / flow.para.dt;
			flow.field.mpg[0] += dmpg1;
			flow.field.mpg[2] += dmpg3;
			for (int j=1; j<flow.para.Ny; j++) {
				flow.field.U.com1.layerAdd(- flow.para.dt * dmpg1, j);
				flow.field.U.com3.layerAdd(- flow.para.dt * dmpg3, j);
			}

# elif CONTROL_TYP == 3
// type 3: dynamically adjust mean pressure gradients to balance the total shear stress calculated from the full-sized channel
			int jb1 = (flow0.para.Ny - flow.para.Ny) / 2, jb2 = flow0.para.Ny - jb1;
			Scla q1(Mesh(flow0.para.Nx, 0, flow0.para.Nz, flow0.para.Lx, 0, flow0.para.Lz)), \
				 q2(q1.meshGet()), q3(q1.meshGet()), q4(q1.meshGet()), \
				 v1(q1.meshGet()), v2(q1.meshGet());
			double dmpg1, dmpg3, tau1, tau3;
			static double mpg1 = 0, mpg3 = 0;

			v1.layerCpy(flow0.field.U.com2, 0, jb1+1).layerAdd(-v1.layerMean());
			v2.layerCpy(flow0.field.U.com2, 0, jb2).layerAdd(-v2.layerMean());


			flow0.field.U.com1.layerUG2CC(q1, 0, jb1);
			flow0.field.U.com1.layerUG2CC(q2, 0, jb1+1);
			flow0.field.U.com1.layerUG2CC(q3, 0, jb2-1);
			flow0.field.U.com1.layerUG2CC(q4, 0, jb2);

			tau1 = .5 / flow0.para.Re * (
				(q2.layerMean() - q1.layerMean()) / flow0.mesh.h[jb1+1] - \
				(q4.layerMean() - q3.layerMean()) / flow0.mesh.h[jb2]	);

			q1.layerMlt(flow0.mesh.dy[jb1+1]);
			q2.layerMlt(flow0.mesh.dy[jb1]);
			q1.layersAdd(q2).layerMlt(.5/flow0.mesh.h[jb1+1]);
			q1.layerAdd(-q1.layerMean());

			q3.layerMlt(flow0.mesh.dy[jb2]);
			q4.layerMlt(flow0.mesh.dy[jb2-1]);
			q3.layersAdd(q4).layerMlt(.5/flow0.mesh.h[jb2]);
			q3.layerAdd(-q3.layerMean());

			tau1 -= .5 * (q1.layersMlt(v1).layerMean() - q3.layersMlt(v2).layerMean());


			flow0.field.U.com3.layerWG2CC(q1, 0, jb1);
			flow0.field.U.com3.layerWG2CC(q2, 0, jb1+1);
			flow0.field.U.com3.layerWG2CC(q3, 0, jb2-1);
			flow0.field.U.com3.layerWG2CC(q4, 0, jb2);

			tau3 = .5 / flow0.para.Re * (
				(q2.layerMean() - q1.layerMean()) / flow0.mesh.h[jb1+1] - \
				(q4.layerMean() - q3.layerMean()) / flow0.mesh.h[jb2]	);

			q1.layerMlt(flow0.mesh.dy[jb1+1]);
			q2.layerMlt(flow0.mesh.dy[jb1]);
			q1.layersAdd(q2).layerMlt(.5/flow0.mesh.h[jb1+1]);
			q1.layerAdd(-q1.layerMean());

			q3.layerMlt(flow0.mesh.dy[jb2]);
			q4.layerMlt(flow0.mesh.dy[jb2-1]);
			q3.layersAdd(q4).layerMlt(.5/flow0.mesh.h[jb2]);
			q3.layerAdd(-q3.layerMean());

			tau3 -= .5 * (q1.layersMlt(v1).layerMean() - q3.layersMlt(v2).layerMean());


			dmpg1 = -2./flow.para.Ly * tau1 - flow.field.mpg[0];
			dmpg3 = -2./flow.para.Ly * tau3 - flow.field.mpg[2];
			flow.field.mpg[0] += dmpg1;
			flow.field.mpg[2] += dmpg3;
			for (int j=1; j<flow.para.Ny; j++) {
				flow.field.U.com1.layerAdd(- flow.para.dt * (dmpg1 - .5*(flow.field.mpg[0]-mpg1)), j);
				flow.field.U.com3.layerAdd(- flow.para.dt * (dmpg3 - .5*(flow.field.mpg[2]-mpg3)), j);
			}
			mpg1 = flow.field.mpg[0];
			mpg3 = flow.field.mpg[2];
		
			q1.meshGet().freeall();

# elif CONTROL_TYP == 4
// type 4: keep the pressure gradients equal to that of the full-sized channel
			double dmpg1 = flow0.field.mpg[0] - flow.field.mpg[0];
			double dmpg3 = flow0.field.mpg[2] - flow.field.mpg[2];
			flow.field.mpg[0] += dmpg1;
			flow.field.mpg[2] += dmpg3;
			for (int j=1; j<flow.para.Ny; j++) {
				flow.field.U.com1.layerAdd(- flow.para.dt * dmpg1, j);
				flow.field.U.com3.layerAdd(- flow.para.dt * dmpg3, j);
			}
# endif

		}


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





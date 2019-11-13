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


double initiate(Field &field, Para &para);
void output(int tstep, double time, Para &para, Field &field, Statis &stas);
void input (int tstep, Para &para, char workpath[1024]);


int main()
{
	Para para("");
	para.showPara();

	Mesh mesh(para.Nx, para.Ny, para.Nz, para.Lx, para.Ly, para.Lz);
	mesh.initYmesh(mesh.getYmesh(para.dy_min));
	mesh.checkYmesh(para.statpath);

	Field field(mesh);
	Statis stas(mesh);

	double time = initiate(field, para);
	int tstep = (time < 1e-10 ? 0 : para.nread); // continue last case or start a new case depending on whether a continue time is read
	
	if (tstep == 0)	output(tstep, time, para, field, stas);

	// main loop
	while (tstep++ < para.Nt) {
		time += para.dt;

		field.bcond(tstep);					// set boundary conditions
		field.getnu(para.Re, para.bftype);	// get the viscosity field
		field.getup(para.dt, para.nthrds);	// time evolution

		if (para.bftype == 1) field.removeSpanMean();	// for MFU
		
		output(tstep, time, para, field, stas);
		input (tstep, para, "");
	}

}




double initiate(Field &field, Para &para)
{
	double time = 0;

	if (para.nread == 0) {
		field.initField(para.inener);
		cout << endl << "Flow fields initiated from laminar." << endl;
	}
	else {
		Para para0(para.inpath);
		cout << endl << "Parameters for continue computing:" << endl;
		para0.showPara();

		Mesh mesh0(para0.Nx, para0.Ny, para0.Nz, para0.Lx, para0.Ly, para0.Lz);
		mesh0.initYmesh(mesh0.getYmesh(para0.statpath));

		Field field0(mesh0);
		field0.readField(para0.fieldpath, para.nread);

		field.initField(field0);
		cout << endl << "Flow fields initiated from existing fields." << endl;

		// read continue time for inplace continuing computation
		if (! strcmp(para0.statpath, para.statpath)) // not the best way to determine whether two paths are quivalent
			time = Statis(mesh0).getLogtime(para0.statpath, para.nread);

		mesh0.freeall();
	}
	return time;
}

void output(int tstep, double time, Para &para, Field &field, Statis &stas)
{
	if (tstep % para.nwrite == 0) {
		field.writeField(para.fieldpath, tstep);
		cout << "Files successfully written for step " << tstep << endl;
	}
	if (tstep % para.nprint == 0) {
		stas.check(field.U, field.P, field.NU, para.Re, para.dt);
		stas.writeProfile(para.statpath);
		stas.writeLogfile(para.statpath, tstep, time);
	}
}

void input(int tstep, Para &para, char workpath[1024])
{
	static time_t last_time = 0;

	if (tstep % para.nprint == 0) {
		char filename[1024]; strcat(strcpy(filename, workpath), "XINDAT");
		struct stat buf;
		stat(filename, &buf);

		if ((int)last_time == 0) last_time = buf.st_mtime;

		if (buf.st_mtime > last_time) {
			last_time = buf.st_mtime;

			para.readPara(workpath);
			cout << endl << filename << " updated at step " << tstep << " as:" << endl;
			para.showPara();
		}
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





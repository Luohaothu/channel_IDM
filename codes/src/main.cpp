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
void input (int tstep, Para &para);



// Vctr& BC_generator(double bc_time)
// {
// 	static Para para("bcpath/XINDAT");
// 	static Mesh mesh(para.Nx, para.Ny, para.Nz);
// 	if (mesh.checkYmesh(mesh.y)) mesh.initYmesh (para.Lx, para.Ly, para.Lz, mesh.getYmesh(para.dy_min, para.Ly));
// 	static Field field(mesh);
// 	static double time = initiate(field, para);
// 	static int tstep = (time < 1e-10 ? 0 : para.nread);

// 	while (bc_time - time > 1e-10) {
// 		time += para.dt;

// 		field.bcond();
// 		field.getnu(para.Re, 0);
// 		field.getup(para.dt, para.nthrds);
// 	}

// 	return field.UBC;
// }

int main()
{
	Para para("XINDAT");
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
		input (tstep, para);
	}

}




double initiate(Field &field, Para &para)
{
	if (para.nread == 0) {
		field.initField(para.inener);
		cout << endl << "Flow fields initiated from laminar." << endl;
		return 0;
	}

	char str[1024];

	sprintf(str, "%s%s", para.inpath, "XINDAT");
	Para para0(str);
	cout << endl << "Parameters for continue computing:" << endl;
	para0.showPara();

	sprintf(str, "%s%s", para.inpath, para0.statpath);
	Mesh mesh0(para0.Nx, para0.Ny, para0.Nz, para0.Lx, para0.Ly, para0.Lz);
	mesh0.initYmesh(mesh0.getYmesh(str));

	sprintf(str, "%s%s", para.inpath, para0.fieldpath);
	Field field0(mesh0);
	field0.readField(str, para.nread);
	field.initField(field0);
	cout << endl << "Flow fields initiated from existing fields." << endl;

	sprintf(str, "%s%s", para.inpath, para0.statpath);
	Statis stas0(mesh0);
	if (! strcmp(para.inpath, "")) // only read continue time for inplace (inpath="") continue computation
		return stas0.getLogtime(str, para.nread);
	else return 0;
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

void input(int tstep, Para &para)
{
	static time_t last_time = 0;
	char filename[] = "XINDAT";

	if (tstep % para.nprint == 0) {
		struct stat buf;
		stat(filename, &buf);

		if ((int)last_time == 0) last_time = buf.st_mtime;

		if (buf.st_mtime > last_time) {
			last_time = buf.st_mtime;

			para.readPara(filename);
			cout << endl << "Parameters updated at step " << tstep << " as:" << endl;
			para.showPara();
		}
	}
}












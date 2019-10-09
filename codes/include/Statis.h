# pragma once

# include "Mesh.h"
# include "Field.h"


class Statis
{
	public:
		Statis(int dim[3]);
		~Statis();

		void checkStat(double *UP[4], class Mesh *pmesh, class Field *pfield, double Re, double dt)
		{
			double *u = UP[0], *v = UP[1], *w = UP[2], *p = UP[3];
			this->checkDiv(u, v, w, pmesh);
			this->checkCFL(u, v, w, pmesh, dt);
			this->checkMean(u, v, w, p, pmesh);
			this->checkTauw(pmesh, Re);
			this->checkEner(u, v, w, p, pmesh, pfield);
		};
		
		double checkDiv(double *u, double *v, double *w, class Mesh *pmesh);
		double checkCFL(double *u, double *v, double *w, class Mesh *pmesh, double dt);
		double checkMean(double *u, double *v, double *w, double *p, class Mesh *pmesh);
		double checkTauw(class Mesh *pmesh, double Re);
		double checkEner(double *u, double *v, double *w, double *p, class Mesh *pmesh, class Field *pfield);

		void writeProfile(char *path, int tstep, double *yc);
		void writeLogfile(char *path, int tstep, double time);
		long int readLogfile(char *path, int tstep, double *time = NULL);

	private:
		const int Nx, Ny, Nz, Nxz;

		double div;	int divpos[3];
		double cfl;	int cflpos[3];
		double velm[3];
		double tauw[3];
		double ener;

		double *Um, *Vm, *Wm, *Pm;
		double *R11,*R22,*R33,*R12,*R23,*R13;
		double *Rpu,*Rpv,*Rpw,*Rpp;
};
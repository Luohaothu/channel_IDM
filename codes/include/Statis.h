# pragma once

# include "Mesh.h"
# include "Field.h"


class Statis
{
	public:
		Statis(int *dim);
		~Statis();

		void checkStat(double **UP, class Mesh *pmesh, class Field *pfield, double Re, double dt)
		{
			this->checkDiv(UP, pmesh);
			this->checkCFL(UP, pmesh, dt);
			this->checkMean(UP, pmesh);
			this->checkTauw(pmesh, Re);
			this->checkEner(UP, pmesh, pfield);
		};
		
		double checkDiv(double **U, class Mesh *pmesh);
		double checkCFL(double **U, class Mesh *pmesh, double dt);
		double checkMean(double **UP, class Mesh *pmesh);
		double checkTauw(class Mesh *pmesh, double Re);
		double checkEner(double **UP, class Mesh *pmesh, class Field *pfield);

		void writeProfile(char *path, int tstep, double *yc);
		void writeLogfile(char *path, int tstep, double time);
		long int readLogfile(char *path, int tstep, double *time = NULL);

	private:
		int Nx, Ny, Nz, Nxz;

		double div;	int divpos[3];
		double cfl;	int cflpos[3];
		double velm[3];
		double tauw[3];
		double ener;

		double *Um, *Vm, *Wm, *Pm;
		double *R11,*R22,*R33,*R12,*R23,*R13;
		double *Rpu,*Rpv,*Rpw,*Rpp;
};
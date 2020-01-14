# pragma once

# include "Basic.h"

class Statis: private Mesh
{
	public:
		Statis(const Mesh &mesh);
		~Statis();

		void check(Vctr &U, Scla &P, Scla &NU, double Re, double dt)
		{
			this->checkDiv (U);
			this->checkCFL (U, dt);
			this->checkTaub(U, Re);
			this->checkEner(U, P, NU);
		};
		double checkDiv (Vctr &U);
		double checkCFL (Vctr &U, double dt);
		double checkTaub(Vctr &U, double Re);
		double checkEner(Vctr &U, Scla &P, Scla &NU);
		
		void writeProfile(char *path, int tstep=-1);

		long int getLogpos(char *path, int tstep);
		double getLogtime(char *path, int tstep);
		void writeLogfile(char *path, int tstep, double time, double mpg[3]);

	private:
		double div;	int divpos[3];
		double cfl;	int cflpos[3];
		double velm[3];	// bulk mean velocities in 3 directions
		double taub[3];	// 3 components of total stress acting on the boundary
		double ener;	// total fluctuation energy in the domain

		double *Um,  *Vm,  *Wm,  *Pm;
		double *R11, *R22, *R33, *R12, *R23, *R13;
		double *Rpu, *Rpv, *Rpw, *Rpp;
		double *Num;
};
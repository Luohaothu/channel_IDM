# pragma once

# include "Basic.h"

class Statis: private Mesh
{
	public:
		Statis(const Mesh &mesh);
		~Statis();

		void check(Vctr &U, Scla &P, Scla &NU, double Re, double dt)
		{
			this->checkDiv(U);
			this->checkCFL(U, dt);
			this->checkMean(U, P, NU);
			this->checkTauw(Re);
			this->checkEner(U, P);
		};
		double checkDiv(Vctr &U);
		double checkCFL(Vctr &U, double dt);
		double checkMean(Vctr &U, Scla &P, Scla &NU);
		double checkTauw(double Re);
		double checkEner(Vctr &U, Scla &P);
		
		void writeProfile(char *path, int tstep=-1);

		long int getLogpos(char *path, int tstep);
		double getLogtime(char *path, int tstep);
		void writeLogfile(char *path, int tstep, double time);

	private:
		double div;	int divpos[3];
		double cfl;	int cflpos[3];
		double velm[3];
		double tauw[3];
		double ener;

		double *Um,  *Vm,  *Wm,  *Pm;
		double *R11, *R22, *R33, *R12, *R23, *R13;
		double *Rpu, *Rpv, *Rpw, *Rpp;
		double *Num;
};
# pragma once

# include "Basic.h"

class Statis: private Mesh
{
	public:
		Statis(): Mesh(1,0,1,1,0,1) { freeall(); }; // default constructor, for static function call
		Statis(const Mesh &mesh);
		~Statis();

		void check(const Feld &FLD, const Scla &NU, double Re, double dt)
		{
			checkDiv (FLD.V);
			checkCFL (FLD.V, dt);
			checkTaub(FLD.V, Re);
			checkEner(FLD.V, FLD.S, NU);
		};
		
		void writeProfile(const char *path, int tstep=-1) const;
		void writeLogfile(const char *path, int tstep, double time, const double mpg[3]) const;
		static double getLogtime(const char *path, int tstep);
		
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

		double checkDiv (const Vctr &U);
		double checkCFL (const Vctr &U, double dt);
		double checkTaub(const Vctr &U, double Re);
		double checkEner(const Vctr &U, const Scla &P, const Scla &NU);

		static long int getLogpos(const char *path, int tstep);
};
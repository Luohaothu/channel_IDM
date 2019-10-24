# pragma once

# include "Basic.h"

class IDM: private Mesh
{
	public:
		IDM(const Mesh &mesh): Mesh(mesh) {};

		// configuration
		double dtset(double dt) { return (this->dt = dt); };
		int ompset(int n);

		// computation interfaces
		void uhcalc(Vctr &UH, Vctr &U , Scla &P, Vctr &UBC, Scla &NU, double mpg[3]);
		void dpcalc(Scla &DP, Vctr &UH, Scla &P, Vctr &UBC);
		void upcalc(Vctr &U , Scla &P , Vctr &UH, Scla &DP, double *mpg);

		// computations functions
		// step 1
		void urhs1(double *ruh, double *u, double *v, double *w, double *p, double *nu, double mpg1);
		void urhs2(double *rvh, double *u, double *v, double *w, double *p, double *nu             );
		void urhs3(double *rwh, double *u, double *v, double *w, double *p, double *nu, double mpg3);
		void mbc(
			double *ruh,double *rvh,double *rwh,
			double *u,	double *v,	double *w,
			double *ubc,double *vbc,double *wbc, double *nu	);
		void getuh1(double *uh,							double *u, double *v, double *w, double *nu);
		void getuh2(double *uh, double *vh,				double *u, double *v, double *w, double *nu);
		void getuh3(double *uh, double *vh, double *wh,	double *u, double *v, double *w, double *nu);
		// step 2
		void rhsdp(double *rdp, double *uh, double *vh, double *wh, double *vbc);
		void getfdp(double *fdp, double refp);
		// step 3
		void update(
			double *u,	double *v,	double *w,	double *p,
			double *uh,	double *vh,	double *wh,	double *dp);
		void meanpg(double *mpg, Vctr &U);

	private:
		double dt;
		int nthrds;	// number of threads for openmp
};






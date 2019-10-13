# pragma once

# include "Matrix.h"
# include "Mesh.h"
# include "Field.h"

class IDM
{
	public:
		void initIDM(double Re, double dt, class Mesh *pmesh);
		void ruhcalc(double *RUH[3], double *U[3], double *P[2], double *UBC[3]);
		void uhcalc(double *UH[3], double *U[3]);
		void dpcalc(double *DP[2], double *UH[3], double *UBC[3], class Field *pfield);
		void upcalc(double *U[3], double *P[2], double *UPH[4], class Mesh *pmesh);

	private:
		int Nx, Ny, Nz, Nxz;
		int IDX(int i, int j, int k) {return Nxz * j + Nx * k + i;};

		// Equation parameters
		double Re, dt;

		// Mesh parameters
		double *dy, *h;				// y interval of the V grid and of the U,W,P grid
		double *hm, *hc, *hp;		// coefficients for wall-normal 2nd-order derivative of U,W
		double *dym, *dyc, *dyp;	// coefficients for wall-normal 2nd-order derivative of V
		double *pmj, *pcj, *ppj;	// coefficients for wall-normal 2nd-order derivative of P in Poisson equation
		double dx, dz;				// grid interval of X and Z
		double dx2, dz2;			// square of grid interval of X and Z

		double *ak1, *ak3;			// coefficients of the Fourier transformed Laplacian operator (related to wavenumbers)

		int *kpa, *kma;				// the forward and backward node indices in periodic Z direction
		int *ipa, *ima;				// the forward and backward node indices in periodic X direction

		// matrix operations initiation
		class Matrix *pmatx;
		class Matrix *pmatyu;
		class Matrix *pmatyv;
		class Matrix *pmatz;


		void urhs1(double *ruh, double *u, double *v, double *w, double *p, double mpg1);
		void urhs2(double *rvh, double *u, double *v, double *w, double *p);
		void urhs3(double *rwh, double *u, double *v, double *w, double *p, double mpg3);
		void mbc(
			double *ruh,double *rvh,double *rwh,
			double *u,	double *v,	double *w,
			double *ubc,double *vbc,double *wbc	);
		void getuh1(double *uh,							double *u, double *v, double *w);
		void getuh2(double *uh, double *vh,				double *u, double *v, double *w);
		void getuh3(double *uh, double *vh, double *wh,	double *u, double *v, double *w);

		void rhsdp(double *rdp, double *uh, double *vh, double *wh, double *vbc);
		void getfdp(double *fdp);

		void update(
			double *u,	double *v,	double *w,	double *p,
			double *uh,	double *vh,	double *wh,	double *dp);
		void meanpg(double *mpg, double *u, double *w, class Mesh *pmesh);
};






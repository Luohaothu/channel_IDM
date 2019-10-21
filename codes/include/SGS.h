# pragma once

# include "Mesh.h"

class SGS
{
	public:
		SGS(int *dim): Nx(dim[0]), Ny(dim[1]), Nz(dim[2]), Nxz(dim[0]*dim[2]) { cs = 0.01; };
		~SGS() {};

		void viscalc(double *nu, double *U[3], double Re, class Mesh *pmesh);

	private:
		int Nx, Ny, Nz, Nxz;
		int IDX(int i, int j, int k) {return Nxz * j + Nx * k + i;};

		double cs; // smargorinsky constant

		void smargorinsky(double *ev, double *u, double *v, double *w, double Re, class Mesh *pmesh);
};
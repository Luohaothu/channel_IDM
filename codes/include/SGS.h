# pragma once

# include "Mesh.h"

class SGS
{
	public:
		SGS(int *dim);
		~SGS();

		void smargorinsky();

	private:
		int Nx, Ny, Nz, Nxz;
		int IDX(int i, int j, int k) {return Nxz * j + Nx * k + i;};


		double *ev, *evi, *evj, *evk;
};
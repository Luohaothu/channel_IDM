# pragma once

# include "Basic.h"


class Interp
{
	public:
		Interp(Scla &s0, Scla &s1);
		~Interp();
		
		void layerPrdLin(int j0, int j1);
		void layerPrdFlt(int j0, int j1);
		void layerTriFlt(int j0, int j1);
		void layerY(int j1, char stgtyp='U');
		void interpolate(char stgtyp);

	private:
		Scla &src, &dst;
		Mesh &mesh0, &mesh1;

		const int Nx0, Ny0, Nz0, Nxz0;
		const int Nx1, Ny1, Nz1, Nxz1;
		const double Lx0, Ly0, Lz0, dx0, dz0;
		const double Lx1, Ly1, Lz1, dx1, dz1;
		const double rscl0, rscl1;

		int *i0l, *k0l;
		int *i0f, *k0f, *nif, *nkf;
		int *j0u, *j0v;
};
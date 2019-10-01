# pragma once




class Mesh
{
	public:
		Mesh(int *dim);
		~Mesh();

		int Nx, Ny, Nz, Nxz;
		double Lx, Ly, Lz;			// domain lengths

		double *y, *yc;				// y coordinates in wall-normal direction (y for V, yc for U,W,P)
		double *dy, *h;				// y interval of the V grid and of the U,W,P grid
		double *hm, *hc, *hp;		// coefficients for wall-normal 2nd-order derivative of U,W
		double *dym, *dyc, *dyp;	// coefficients for wall-normal 2nd-order derivative of V
		double *pmj, *pcj, *ppj;	// coefficients for wall-normal 2nd-order derivative of P in Poisson equation
		double dx, dz;				// grid interval of X and Z
		double dx2, dz2;			// square of grid interval of X and Z
		double vol;					// volume of the whole physical domain

		double *ak1, *ak3;			// coefficients of the Fourier transformed Laplacian operator (related to wavenumbers)

		int *kpa, *kma;				// the forward and backward node indices in periodic Z direction
		int *ipa, *ima;				// the forward and backward node indices in periodic X direction

		void initMesh(double *len, char *path, double dy_min);

		// tool functions
		double divergence(double *u, double *v, double *w, int i, int j, int k);
		double convection(double *u, double *v, double *w, int i, int j, int k);
		double* gradient(double *p, int i, int j, int k);

		double yMeanU(double *src, int i, int k);
		double yMeanV(double *src, int i, int k);

		double layerMean(double *src, int j = 0);
		double* layerCenterU(double *dst, double *src, int j1 = 0, int j0 = 0);	// interpolate j0 layer of src from U grid to cell center, result stored in j1 layer of dst
		double* layerCenterV(double *dst, double *src, int j1, int j0);			// interpolate src from V grid to the j0 layer of P grid, result stored in J1 layer of dst
		double* layerCenterW(double *dst, double *src, int j1 = 0, int j0 = 0);	// interpolate j0 layer of src from W grid to cell center, result stored in j1 layer of dst

		double bulkMeanU(double *src);
		double bulkMeanV(double *src);

	private:
		int IDX(int i, int j, int k) {return Nxz * j + Nx * k + i;};
		double getYmesh(char *path, double dy_min);
		double distrib_Y(int j, double gamma);
};





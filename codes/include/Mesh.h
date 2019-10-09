# pragma once


# define PI 3.1415926535898


class Mesh
{
	public:
		Mesh(int dim[3]);
		~Mesh();

		const int Nx, Ny, Nz, Nxz;
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

		void initMesh(double len[3], char *path, double dy_min);
		bool checkYmesh(char *path, double *ymesh);

		double kx(int i) { return ( i - ( i > (int)(Nx/2) ? Nx : 0 ) ) * (2.0*PI/Lx); };
		double kz(int k) { return ( k - ( k > (int)(Nz/2) ? Nz : 0 ) ) * (2.0*PI/Lz); };

		// tool functions
		// vector & tensor operators
		double divergence(double *u, double *v, double *w, int i, int j, int k);
		double convection(double *u, double *v, double *w, int i, int j, int k);
		double* gradient(double *p, int i, int j, int k);
		double* strainrate(double *u, double *v, double *w, int i, int j, int k);
		// integration
		double yMeanU(double *src, int i, int k);
		double yMeanV(double *src, int i, int k);
		double layerMean(double *src, int j = 0);
		double bulkMeanU(double *src);
		double bulkMeanV(double *src);
		// interpolation
		double* layerUG2CC(double *dst, double *src, int j1=0, int j0=0);
		double* layerVG2CC(double *dst, double *src, int j1, int j0);
		double* layerWG2CC(double *dst, double *src, int j1=0, int j0=0);
		double* layerCC2UG(double *dst, double *src, int j1=0, int j0=0);
		double* layerCC2VG(double *dst, double *src, int j1, int j0);
		double* layerCC2WG(double *dst, double *src, int j1=0, int j0=0);


	private:
		int IDX(int i, int j, int k) {return Nxz * j + Nx * k + i;};
		double* getYmesh(char *path);
		double* getYmesh(double dy_min);
		double distrib_Y(int j, double gamma);

};





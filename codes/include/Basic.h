# pragma once

# include <fftw3.h>

# define PI 3.1415926535898


class Bulk
{
	public:
		Bulk(int n1, int n2, int n3, bool inift=false);	// memory allocation, dim is the array dimension, unrelated to the mesh dimension
		void freeall();

		double* fft ();
		double* ifft();

		// memory access functions
		double& id (int i, int j=0, int k=0) { return  q [nxz * j + nx * k + i]; };
		double& idf(int i, int j=0, int k=0) { return fq[nxzr* j + nxr* k + i]; };

		// tool functions
		double* layerGet (int j=0);
		double* layerCpy (double *src, int j1=0, int j0=0);	// copy layer j0 of src to layer j1 of dst
		double* layerSet (double a, int j=0);				// set the value of layer j of dst to a
		double* layerAdd (double a, int j=0);				// add value a to layer j of dst
		double* layerMlt (double a, int j=0);				// multiply value a to layer j of dst
		double* layersAdd(double *src, int j1=0, int j0=0);	// add layer j0 of src to layer j1 of dst
		double* layersMlt(double *src, int j1=0, int j0=0);	// multiply layer j0 of src to layer j1 of dst

		double* bulkGet () { return q; };
		double* bulkGetF() { return fq; };
		double* bulkCpy (double *src);
		double* bulkSet (double a);
		double* bulkAdd (double a);
		double* bulkMlt (double a);

		Bulk& layerCpy (Bulk &b, int j1=0, int j0=0) { this->layerCpy (b.bulkGet(), j1, j0); return *this; };
		Bulk& layersAdd(Bulk &b, int j1=0, int j0=0) { this->layersAdd(b.bulkGet(), j1, j0); return *this; };
		Bulk& layersMlt(Bulk &b, int j1=0, int j0=0) { this->layersMlt(b.bulkGet(), j1, j0); return *this; };
		Bulk& bulkCpy  (Bulk &b) { this->bulkCpy(b.bulkGet()); return *this; };

		// IO functions
		void fileIO(char *path, char *name, char mode);
		void debug_AsciiOutput(char *path, char *name, int j1, int j2);

	protected:
		double *q;
		double *fq;		// Fourier tansformed q with only the kx >= 0 half

	private:
		const int nx, ny, nz, nxz;
		const int nxc, nxr, nxzc, nxzr;
		fftw_plan *frcs;	// list of fft plans: 2D FFT from real data to complex data, FORWARD (exponent -1) and no normalization
		fftw_plan *fcrs;	// list of ifft plans: 2D FFT from complex data to real data, BACKWARD (exponent 1) with normalization
};


class Mesh
{
	public:
		Mesh(int n1, int n2, int n3);	// memory allocation, dim is the mesh dimension, whose Ny is 1 less than the array dimension
		void freeall();

		const int Nx, Ny, Nz, Nxz;
		int IDX(int i, int j, int k) { return Nxz * j + Nx * k + i; };

		double Lx, Ly, Lz;			// domain lengths
		double dx, dz, dx2, dz2;	// grid interval of X & Z and their squares
		double vol, *dvol;			// volume of the whole physical domain and of every cell

		double *y,  *yc;			// y coordinates in wall-normal direction (y for V, yc for U,W,P)
		double *dy, *h;				// y interval of the wall-normal grid (dy for V, h for U,W,P)
		double *hm,  *hc,  *hp;		// coefficients for wall-normal 2nd-order derivative of U,W
		double *dym, *dyc, *dyp;	// coefficients for wall-normal 2nd-order derivative of V

		double *pmj, *pcj, *ppj;	// coefficients for wall-normal 2nd-order derivative of P in Poisson equation, boundaries excluded
		double *ak1, *ak3;			// coefficients of the Fourier transformed Laplacian operator (related to wavenumbers)

		int *kpa, *kma;				// the forward and backward node indices in periodic Z direction
		int *ipa, *ima;				// the forward and backward node indices in periodic X direction

		double kx(int i) { return ( i - ( i > (int)(Nx/2) ? Nx : 0 ) ) * (2.0*PI/Lx); };
		double kz(int k) { return ( k - ( k > (int)(Nz/2) ? Nz : 0 ) ) * (2.0*PI/Lz); };

		void initMesh(double l1, double l2, double l3, double *ymesh);
		void checkYmesh(double *ymesh, char *path=NULL);
		double* getYmesh(double dy_min, double y_max);
		double* getYmesh(char *path);
};


// # include <iostream>
// using namespace std;

class Scla: public Mesh, public Bulk
{
	public:
		Scla(const Mesh &mesh, const Bulk &bulk): Mesh(mesh), Bulk(bulk) {};
		Scla(const Mesh &mesh, bool inift=false): Mesh(mesh), Bulk(Nx, Ny+1, Nz, inift) {};

		// Mesh& meshGet() { return (Mesh&)(*this); };

		// integration
		double yMeanU(int i, int k);
		double yMeanV(int i, int k);
		double layerMean(int j=0);
		double bulkMeanU();
		double bulkMeanV();

		// interpolation
		double* layerUG2CC(double *src, int j1=0, int j0=0);
		double* layerVG2CC(double *src, int j1,   int j0);
		double* layerWG2CC(double *src, int j1=0, int j0=0);

		Scla& layerUG2CC(Scla &s, int j1=0, int j0=0) { this->layerUG2CC(s.bulkGet(), j1, j0); return *this; };
		Scla& layerVG2CC(Scla &s, int j1,   int j0)   { this->layerVG2CC(s.bulkGet(), j1, j0); return *this; };
		Scla& layerWG2CC(Scla &s, int j1=0, int j0=0) { this->layerWG2CC(s.bulkGet(), j1, j0); return *this; };

		Scla& interpolate(Scla &src);

		// operators
		double* gradient(int i, int j, int k);
};


class Vctr: private Mesh
{
	public:
		Vctr(const Mesh &mesh, const Scla &s1, const Scla &s2, const Scla &s3):
			Mesh(mesh),
			com1(s1), com2(s2), com3(s3),
			u(com1.bulkGet()), v(com2.bulkGet()), w(com3.bulkGet()) {};

		Vctr(const Mesh &mesh, bool inift=false):
			Mesh(mesh),
			com1(mesh, inift), com2(mesh, inift), com3(mesh, inift),
			u(com1.bulkGet()), v(com2.bulkGet()), w(com3.bulkGet()) {};

		Scla com1, com2, com3;

		double* bulkGet(int i) { return (i==1 ? u : i==2 ? v : i==3 ? w : NULL); };
		
		// void interpolate(Vctr &src) {
		// 	com1.interpolate(src.com1);
		// 	com2.interpolate(src.com2);
		// 	com3.interpolate(src.com3);
		// };

		// vector operators
		double  divergence(int i, int j, int k);
		double  convection(int i, int j, int k);
		double* strainrate(int i, int j, int k);

	private:
		double *u;
		double *v;
		double *w;
};




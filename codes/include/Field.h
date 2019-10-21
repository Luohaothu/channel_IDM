# pragma once

# include <string>
# include <math.h>
# include <fftw3.h>
typedef fftw_complex fcmplx;

# include "Mesh.h"


class Field
{
	public:
		Field(int dim[3]);	// allocate memory for pointers
		~Field();

		double *U[3], *P[2], *UP[4];
		double *UH[3], *DP[2], *UPH[4];
		double *UBC[3];
		double *nu;	// viscosity field

		void initField(double energy, class Mesh *pmesh);	// initiate flow fields from laminar with random fluctuations
		void initField(class Field *pf0, class Mesh *pm0, class Mesh *pm);	// initiate flow fields from existing fields
		void bcond(int tstep);

		void applyBC();
		void applyBC(double dt);

		// void bodyForce(int bftype);


		void fft()	{ for (int j=0; j<=Ny; j++) fftw_execute(frcs[j]); };
		void ifft()	{ for (int j=0; j<=Ny; j++) fftw_execute(fcrs[j]);	this->bulkMult(dp, 1.0/Nxz); };

		// memory access functions
		int IDX(int i, int j, int k) {return Nxz * j + Nx * k + i;};
		int IDXF(int i, int j, int k) {return Nxzr * j + Nxr * k + 2*i;};

		// tool functions
		double* layerCopy	(double *dst, double *src, int j1 = 0, int j0 = 0);	// copy layer j0 of src to layer j1 of dst
		double* layerSet	(double *dst, double a, int j = 0);	// set the value of layer j of dst to a
		double* layerAdd	(double *dst, double a, int j = 0);	// add value a to layer j of dst
		double* layerMult	(double *dst, double a, int j = 0);	// multiply value a to layer j of dst
		double* layersAdd	(double *dst, double *src, int j1 = 0, int j0 = 0);	// add layer j0 of src to layer j1 of dst
		double* layersMult	(double *dst, double *src = NULL, int j1 = 0, int j0 = 0);	// multiply layer j0 of src to layer j1 of dst. NULL src will be set to dst, which needs caution

		double* bulkCopy	(double *dst, double *src);
		double* bulkSet		(double *dst, double a);
		double* bulkAdd		(double *dst, double a);
		double* bulkMult	(double *dst, double a);

		double* removeSpanMean(double *dst, int j);

		// IO functions
		void writeField(char *path, int tstep);
		void writeFieldDt(char *path, int tstep); // note: the output may not be correct at step 0 
		void readField (char *path, int tstep);
		void fieldIO(char *path, double *ptr, char *name, char mode);
		void writeTecplot(char *path, int tstep, double time, class Mesh *pmesh);
		
		// debug functions
		void debug_AsciiOutput(int tstep, double *ptr, int j1, int j2, std::string filepre);
		void debug_AsciiOutputF(int tstep, int j1, int j2);
		void debug_Output(int tstep);

	private:
		const int Nx, Ny, Nz, Nxz;
		const int Nxc, Nxr, Nxzc, Nxzr;
		fftw_plan *frcs;	// list of fft plans: 2D FFT from real data to complex data, FORWARD (exponent -1) and no normalization
		fftw_plan *fcrs;	// list of ifft plans: 2D FFT from complex data to real data, BACKWARD (exponent 1) and no normalization

		double *u, *v, *w, *p, mpg[3];	// velocity and pressure fields
		double *uh, *vh, *wh, *dp;		// intermediate velocity fields and incremental pressure field
		double *ubc, *vbc, *wbc;		// Dirichlet velocity boundary conditions
		double *fdp;					// Fourier tansformed DP

		/* note:
			uh,vh,wh also serve to store the RHS of delta u^* (increment of intermediate velocities), and delta u^* itself
			dp also serve to store the RHS of the pressure Poisson equation (in physical space)
			fdp also serve to store the RHS of the pressure Poisson equation (in Fourier space)	*/

};





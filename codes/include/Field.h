# pragma once

# include <string>
# include <math.h>
# include <fftw3.h>
typedef fftw_complex fcmplx;

# include "Mesh.h"


class Field
{
	public:
		Field(int *dim);	// allocate memory for pointers
		~Field();

		void initField(double energy, class Mesh *pmesh);
		void bcond(int tstep);
		void applyBC();

		void fft()	{ for (int j=0; j<Ny-1; j++) fftw_execute(frcs[j]); };
		void ifft()	{ for (int j=0; j<Ny-1; j++) fftw_execute(fcrs[j]);	this->bulkMult(dp, 1.0/Nxz); };

		// memory access functions
		double** U	() { static double *U	[3] = {	u , v , w	};	return (double**) U;	};
		double** UP	() { static double *UP	[4] = {	u, v, w, p	};	return (double**) UP;	};
		double** UH	() { static double *UH	[3] = {	uh, vh, wh	};	return (double**) UH;	};
		double** UBC() { static double *UBC	[3] = {	ubc,vbc,wbc	};	return (double**) UBC;	};
		double** P	() { static double *P	[2] = {	p	,	mpg	};	return (double**) P;	};
		double** DP	() { static double *DP	[2] = {	dp	,	fdp	};	return (double**) DP;	};

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

		// IO functions
		void writeField(char *path, int tstep);
		void readField (char *path, int tstep);
		void fieldIO(char *path, double *ptr, char *name, char mode);
		void writeTecplot(char *path, int tstep, double time, class Mesh *pmesh);
		
		// debug functions
		void debug_AsciiOutput(int tstep, double *ptr, int j1, int j2, std::string filepre);
		void debug_AsciiOutputF(int tstep, int j1, int j2);
		void debug_Output(int tstep);

	private:
		int Nx, Ny, Nz, Nxz;
		int IDX(int i, int j, int k) {return Nxz * j + Nx * k + i;};

		int Nxc, Nxr, Nxzc, Nxzr;
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





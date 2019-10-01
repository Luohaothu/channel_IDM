# include <iostream>
# include <stdio.h>
# include <math.h>


# include "Mesh.h"

using namespace std;



Mesh::Mesh(int *dim):
/* allocate memory for pointers */
Nx(dim[0]), Ny(dim[1]), Nz(dim[2]), Nxz(dim[0]*dim[2])
{
	y = new double [Ny+1];	// Ny layers (including 2 walls) + 1 redundant layer (index 0)
	yc= new double [Ny+1];	// Ny+1 layers (including 2 walls)

	dy= new double [Ny+1];	// Ny-1 boxes + 2 redundant (index 0 and Ny)
	h = new double [Ny+1];	// Ny layers + 1 redundant layer (index 0), aligned to y

	hm = new double [Ny];
	hc = new double [Ny];
	hp = new double [Ny];

	dym = new double [Ny];
	dyc = new double [Ny];
	dyp = new double [Ny];

	pmj = new double [Ny];
	pcj = new double [Ny];
	ppj = new double [Ny];

	ak1 = new double [Nx];
	ak3 = new double [Nz];

	kpa = new int [Nz];
	kma = new int [Nz];

	ipa = new int [Nx];
	ima = new int [Nx];
}

Mesh::~Mesh()
{
	delete [] y; delete [] yc;
	delete [] dy; delete [] h; 
	delete [] hm; delete [] hc; delete [] hp;
	delete [] dym; delete [] dyc; delete [] dyp;
	delete [] pmj; delete [] pcj; delete [] ppj;
	delete [] ak1; delete [] ak3;
	delete [] kpa; delete [] kma;
	delete [] ipa; delete [] ima;
}


void Mesh::initMesh(double *len, char *path, double dy_min)
{
	int i, j, k;

	Lx = len[0];
	Ly = len[1];
	Lz = len[2];
	dx = (double) Lx / Nx;
	dz = (double) Lz / Nz;
	dx2 = dx * dx;
	dz2 = dz * dz;
	vol = Lx * Ly * Lz;
	
	// y grids
	if ( abs( this->getYmesh(path, dy_min) - Ly ) > 1e-10 ) printf("Y grid error !!!");

	yc[0] = 0.0;
	yc[Ny] = Ly;
	for (j=1; j<Ny; j++)	yc[j] = 0.5 * ( y[j] + y[j+1] );

	// y intervals
	for (j=1; j<Ny; j++)	{ dy[j] = y[j+1] - y[j]; }
	for (j=2; j<Ny; j++)	{ h[j] = (y[j+1] - y[j-1]) / 2.0; }
	h[1] = dy[1] / 2.0;
	h[Ny] = dy[Ny-1] / 2.0;
	h[0] = 0.0;
	dy[0] = 0.0;
	dy[Ny] = 0.0;

	// coefficients for wall-normal 2nd-order derivative
	for (j=1; j<Ny; j++) {	// for U & W in laplacian operator, boundaries not excluded
		hp[j] = 1.0 / dy[j] / h[j+1];
		hc[j] = 1.0 / dy[j] * ( 1.0/h[j+1] + 1.0/h[j] );	// need to multiply by -1 in use
		hm[j] = 1.0 / dy[j] / h[j];
	}
	for (j=2; j<Ny; j++) {	// for V in laplacian operator, boundaries not excluded
		dyp[j] = 1.0 / h[j] / dy[j];
		dyc[j] = 1.0 / h[j] * ( 1.0/dy[j] + 1.0/dy[j-1] );	// need to multiply by -1 in use
		dym[j] = 1.0 / h[j] / dy[j-1];
	}
	for (j=1; j<Ny; j++) {	// for P in divergence-gradient operator, boundaries excluded
		if (j>1 && j<Ny-1) {
			ppj[j] = hp[j];
			pcj[j] = - hc[j];								// no need to multiply by -1 in use
			pmj[j] = hm[j];		}
		else if (j==1) {
			ppj[j] = hp[1];
			pcj[j] = - 1.0 / dy[1] / h[2];
			pmj[j] = 0.0;		}
		else if (j==Ny-1) {
			ppj[j] = 0.0;
			pcj[j] = - 1.0 / dy[Ny-1] / h[Ny-1];
			pmj[j] = hm[Ny-1];	}
	}

	// initiate wavenumbers
	const double PI = 3.1415926535898;
	for (i=0; i<Nx; i++) ak1[i] = 2.0/dx2 * ( 1.0 - cos( 2*PI/Nx * (i - (i>Nx/2?Nx:0)) ) );
	for (k=0; k<Nz; k++) ak3[k] = 2.0/dz2 * ( 1.0 - cos( 2*PI/Nz * (k - (k>Nz/2?Nz:0)) ) );

	// initiate indexes
	for (k=0; k<Nz; k++) {
		kpa[k] = (k+1) % Nz;
		kma[k] = (k-1+Nz) % Nz;
	}
	for (i=0; i<Nx; i++) {
		ipa[i] = (i+1) % Nx;
		ima[i] = (i-1+Nx) % Nx;
	}
}


double Mesh::getYmesh(char *path, double dy_min)
{
	FILE *fp;
	char str[1024];
	int j;

	sprintf(str, "%sCHANNEL.GRD", path?path:"");

	if (dy_min <= 0) {
	/* read grid from file at specified path */
		fp = fopen(str, "r");
			for (j=0; j<=Ny; j++) {
				fgets(str, 1024, fp);
				sscanf(str, "%le", & y[j]);
			}
		fclose(fp);
	}
	else {
	/* generate hyperbolic tangent grid (Newton iteration to determine parameter) */
		double gamma = 1.0, dgamma = 0.1;
		double F, F0 = this->distrib_Y(2, gamma-dgamma) - dy_min;
		double grad;

		for (j=0; j<100; j++) {
			F = this->distrib_Y(2, gamma) - dy_min; // target equation: F(gamma) = distrib_Y(2,gamma) - dy_min = 0
			grad = (F-F0) / dgamma;	if (! grad) {break;}
			dgamma = - F / grad;	if (! dgamma) {break;}
			gamma += dgamma;
			F0 = F;
		}

		y[0] = 0.0;
		for (j=1; j<=Ny; j++) y[j] = this->distrib_Y(j, gamma);

		// write grid file to specified path
		fp = fopen(str, "w");
			for (j=0; j<=Ny; j++) {
				sprintf(str, "%.18e\n", y[j]);
				fputs(str, fp);
			}
		fclose(fp);
	}
	return y[Ny];
}


double Mesh::distrib_Y(int j, double gamma)
/* generate y grids with hyperbolic tangent distribution */
{
	double delta = Ly / 2.0;
	double ytild = 2.0 * (j-1) / (Ny-1) - 1.0;
	return delta * ( 1.0 + tanh(gamma * delta * ytild) / tanh(gamma * delta) );
}




/********** tool functions **********/

double* Mesh::layerCenterU(double *dst, double *src, int j1, int j0)
{
	int i, k;
	for (k=0; k<Nz; k++) {	// note: pointer dst must NOT be the same as pointer src
	for (i=0; i<Nx; i++) {	dst[IDX(i,j1,k)] = 0.5 * ( src[IDX(i,j0,k)] + src[IDX(ipa[i],j0,k)] );	}}
	return dst;
}
double* Mesh::layerCenterV(double *dst, double *src, int j1, int j0)
{
	int i, k, jd = (j0==0 ? 1 : j0), ju = (j0==Ny ? Ny : j0+1);
	for (k=0; k<Nz; k++) {	// note: pointer dst must NOT be the same as pointer src
	for (i=0; i<Nx; i++) {	dst[IDX(i,j1,k)] = 0.5 * ( src[IDX(i,jd,k)] + src[IDX(i,ju,k)] );	}}
	return dst;
}
double* Mesh::layerCenterW(double *dst, double *src, int j1, int j0)
{
	int i, k;
	for (k=0; k<Nz; k++) {	// note: pointer dst must NOT be the same as pointer src
	for (i=0; i<Nx; i++) {	dst[IDX(i,j1,k)] = 0.5 * ( src[IDX(i,j0,k)] + src[IDX(i,j0,kpa[k])] );	}}
	return dst;
}


double Mesh::layerMean(double *src, int j)
{
	int k, i, kk, jj = Nxz * j;
	double integral = 0.0, dA = dx * dz, A = Lx * Lz;
	for (k=0; k<Nz; k++) { kk = jj + Nx * k;
	for (i=0; i<Nx; i++) { integral += src[kk + i] * dA;	}}
	return integral / A;
}
double Mesh::bulkMeanU(double *src) // only for U, W
{
	double integral = 0.0;
	for (int j=1; j<Ny; j++)	integral += this->layerMean(src, j) * dy[j];
	return integral / Ly;
}
double Mesh::bulkMeanV(double *src) // only for V
{
	double integral = 0.0;
	for (int j=1; j<=Ny; j++)	integral += this->layerMean(src, j) * h[j];
	return integral / Ly;
}
double Mesh::yMeanU(double *src, int i, int k)
{
	int ki = Nx * k + i;
	double integral = 0.0;
	for (int j=1; j<Ny; j++)	integral += src[Nxz * j + ki] * dy[j];
	return integral / Ly;
}
double Mesh::yMeanV(double *src, int i, int k)
{
	int ki = Nx * k + i;
	double integral = 0.0;
	for (int j=1; j<=Ny; j++)	integral += src[Nxz * j + ki] * h[j];
	return integral / Ly;
}


double Mesh::divergence(double *u, double *v, double *w, int i, int j, int k)
/* compute divergence of a vector field at the center of cell (i,j,k) */
{
	int idx= IDX(i,j,k);
	int ip = IDX(ipa[i],j,k);
	int jp = IDX(i,j+1,k);	// 1 <= j <= Ny-1
	int kp = IDX(i,j,kpa[k]);
	return	( u[ip] - u[idx] ) / dx
		+	( v[jp] - v[idx] ) / dy[j]
		+	( w[kp] - w[idx] ) / dz;
}

double Mesh::convection(double *u, double *v, double *w, int i, int j, int k)
/* compute convection coefficient of a vector field at the center of cell (i,j,k) */
{
	int idx= IDX(i,j,k);
	int ip = IDX(ipa[i],j,k);
	int jp = IDX(i,j+1,k);	// 1 <= j <= Ny-1
	int kp = IDX(i,j,kpa[k]);
	return	0.5 * ( u[ip] + u[idx] ) / dx
		+	0.5 * ( v[jp] + v[idx] ) / dy[j]
		+	0.5 * ( w[kp] + w[idx] ) / dz;
}

double* Mesh::gradient(double *p, int i, int j, int k)
/* compute gradient of a scalar field at grid points corresponding to U,V,W */
{
	int idx= IDX(i,j,k);
	int im = IDX(ima[i],j,k);
	int jm = IDX(i,j-1,k);
	int km = IDX(i,j,kma[k]);
	double *grad = new double [3];
	grad[0] = ( p[idx] - p[im] ) / dx;	// 1 <= j <= Ny-1
	grad[1] = ( p[idx] - p[jm] ) / h[j];// 2 <= j <= Ny-1 (j==1 may not be valid)
	grad[2] = ( p[idx] - p[km] ) / dz;	// 1 <= j <= Ny-1
	return grad;
}









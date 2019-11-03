# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>


# include "Basic.h"

using namespace std;


/********** integration functions **********/

double Scla::yMeanU(int i, int k)
{
	int j, ki = Nx * k + i;
	double integral = 0.0;
	for (j=1; j<Ny; j++) integral += this->q[Nxz * j + ki] * dy[j];
	return integral / Ly;
}
double Scla::yMeanV(int i, int k)
{
	int j, ki = Nx * k + i;
	double integral = 0.0;
	for (j=1; j<=Ny; j++) integral += this->q[Nxz * j + ki] * h[j];
	return integral / Ly;
}
double Scla::layerMean(int j)
{
	int k, i, kk, jj = Nxz * j;
	double integral = 0.0, dA = dx * dz, A = Lx * Lz;
	for (k=0; k<Nz; k++) { kk = jj + Nx * k;
	for (i=0; i<Nx; i++) { integral += this->q[kk + i] * dA;	}}
	return integral / A;
}
double Scla::bulkMeanU() // only for U, W
{
	double integral = 0.0;
	for (int j=1; j<Ny; j++) integral += this->layerMean(j) * dy[j];
	return integral / Ly;
}
double Scla::bulkMeanV() // only for V
{
	double integral = 0.0;
	for (int j=1; j<=Ny; j++) integral += this->layerMean(j) * h[j];
	return integral / Ly;
}


/********** interpolation functions **********/

double* Scla::layerUG2CC(double *dst, int j1, int j0)	// CAUTION: pointer dst must NOT be self
/* interpolate j0 layer from U grid of src to j1 layer of cell center of dst */
{
	int i, k, kk0, kk1;
	for (k=0; k<Nz; k++) { kk0 = Nxz * j0 + Nx * k; kk1 = Nxz * j1 + Nx * k;
	for (i=0; i<Nx; i++) { dst[kk1+i] = 0.5 * (q[kk0+i] + q[kk0+ipa[i]]); }}
	return dst;
}
double* Scla::layerVG2CC(double *dst, int j1, int j0)	// CAUTION: pointer dst must NOT be self
/* interpolate j0,j0+1 layers from V grid of src to j1 layer of cell center of dst */
{
	int i, k, kkd, kku, kk1, jd = (j0==0 ? 1 : j0), ju = (j0==Ny ? Ny : j0+1);
	for (k=0; k<Nz; k++) { kkd = Nxz * jd + Nx * k; kku = Nxz * ju + Nx * k; kk1 = Nxz * j1 + Nx * k;
	for (i=0; i<Nx; i++) { dst[kk1+i] = 0.5 * (q[kkd+i] + q[kku+i]);	}}
	return dst;
}
double* Scla::layerWG2CC(double *dst, int j1, int j0)	// CAUTION: pointer dst must NOT be self
/* interpolate j0 layer from W grid of src to j1 layer of cell center of dst */
{
	int i, k, kkb, kkf, kk1;
	for (k=0; k<Nz; k++) { kkb = Nxz * j0 + Nx * k; kkf = Nxz * j0 + Nx * kpa[k]; kk1 = Nxz * j1 + Nx * k;
	for (i=0; i<Nx; i++) { dst[kk1+i] = 0.5 * (q[kkb+i] + q[kkf+i]);	}}
	return dst;
}
// double* Scla::layerUG2CC(double *src, int j1, int j0)	// CAUTION: pointer src must NOT be self
// /* interpolate j0 layer of src from U grid to cell center, result stored in j1 layer of q */
// {
// 	int i, k, kk0, kk1;
// 	for (k=0; k<Nz; k++) { kk0 = Nxz * j0 + Nx * k; kk1 = Nxz * j1 + Nx * k;
// 	for (i=0; i<Nx; i++) { this->q[kk1+i] = 0.5 * (src[kk0+i] + src[kk0+ipa[i]]); }}
// 	return this->q;
// }
// double* Scla::layerVG2CC(double *src, int j1, int j0)	// CAUTION: pointer src must NOT be self
// /* interpolate j0,j0+1 layers of src from V grid to cell center, result stored in j1 layer of q */
// {
// 	int i, k, kkd, kku, kk1, jd = (j0==0 ? 1 : j0), ju = (j0==Ny ? Ny : j0+1);
// 	for (k=0; k<Nz; k++) { kkd = Nxz * jd + Nx * k; kku = Nxz * ju + Nx * k; kk1 = Nxz * j1 + Nx * k;
// 	for (i=0; i<Nx; i++) { this->q[kk1+i] = 0.5 * (src[kkd+i] + src[kku+i]);	}}
// 	return this->q;
// }
// double* Scla::layerWG2CC(double *src, int j1, int j0)	// CAUTION: pointer src must NOT be self
// /* interpolate j0 layer of src from W grid to cell center, result stored in j1 layer of q */
// {
// 	int i, k, kkb, kkf, kk1;
// 	for (k=0; k<Nz; k++) { kkb = Nxz * j0 + Nx * k; kkf = Nxz * j0 + Nx * kpa[k]; kk1 = Nxz * j1 + Nx * k;
// 	for (i=0; i<Nx; i++) { this->q[kk1+i] = 0.5 * (src[kkb+i] + src[kkf+i]);	}}
// 	return this->q;
// }

Scla& Scla::interpolate(Scla &src)
/* initiate flow field from given fields, interpolated to the current grid */
{
	int i, j, k, idx, j0, *i0 = new int [Nx], *k0 = new int [Nz];
	Mesh &mesh0 = src.meshGet();
	int Nx0 = mesh0.Nx, Ny0 = mesh0.Ny, Nz0 = mesh0.Nz, Nxz0 = mesh0.Nxz;
	double *yc0 = mesh0.yc;

	// nearest interpolation for x and z directions in Fourier space
	// find out nearest point for every wavenumber k_z
	for (k=0; k<Nz; k++) {	k0[k] = 0;
	for (int k1=1; k1<Nz0; k1++) {
		if (  fabs(kz(k) - mesh0.kz(k1)   )
			< fabs(kz(k) - mesh0.kz(k0[k]))	)	k0[k] = k1;

	}}
	// find out nearesr point for every wavenumber k_x
	for (i=0; i<Nx; i++) {	i0[i] = 0;
	for (int i1=1; i1<Nx0; i1++) {
		if (  fabs(kx(i) - mesh0.kx(i1)   )
			< fabs(kx(i) - mesh0.kx(i0[i]))	)	i0[i] = i1;
	}}

	src.fft();

	// linear interpolations for y direction
	for (j=0; j<=Ny; j++) {
		// find out the range (j0-1 ~ j0) where j fall
		for (j0=0; j0<=Ny0; j0++) { if (yc0[j0] > yc[j]) break; }
		// if (n == 1) { j0 = j==0 ? 0 : (j==1 ? 2 : j0); }

		for (k=0; k<Nz; k++) {
		for (i=0; i<(int)(Nx/2+1); i++) {
			// values at both ends are extended for out-of-range points
			if (j0 == 0) {
				this->idf(2*i,  j,k) = src.idf(2*i0[i],  j0,k0[k]);	// real part
				this->idf(2*i+1,j,k) = src.idf(2*i0[i]+1,j0,k0[k]);	// imaginary part
			}
			else if (j0 == Ny0+1) {
				this->idf(2*i,  j,k) = src.idf(2*i0[i],  j0-1,k0[k]);
				this->idf(2*i+1,j,k) = src.idf(2*i0[i]+1,j0-1,k0[k]);
			}
			// linear interpolation for in-range points
			else {
				this->idf(2*i,j,k) =
					(yc0[j0] - yc[j]) / ( yc0[j0] - yc0[j0-1] ) * src.idf(2*i0[i],j0-1,k0[k])
				+	(yc[j]-yc0[j0-1]) / ( yc0[j0] - yc0[j0-1] ) * src.idf(2*i0[i],j0,  k0[k]);

				this->idf(2*i+1,j,k) =
					(yc0[j0] - yc[j]) / ( yc0[j0] - yc0[j0-1] ) * src.idf(2*i0[i]+1,j0-1,k0[k])
				+	(yc[j]-yc0[j0-1]) / ( yc0[j0] - yc0[j0-1] ) * src.idf(2*i0[i]+1,j0,  k0[k]);
			}
		}}
	}
	this->ifft();
	this->bulkMlt(1.0 * Nxz / Nxz0);

	delete [] i0; delete [] k0;
	return *this;
}


/********** differentiation operators **********/

double* Scla::gradient(int i, int j, int k)
/* compute gradient of a scalar field at cell center to grid points corresponding to U,V,W */
// CAUTION: avoid successive calling to this function, because the static return variable will be overwritten every time
{
	int idx= IDX(i,j,k);
	int im = IDX(ima[i],j,k);
	int jm = IDX(i,j-1,k);
	int km = IDX(i,j,kma[k]);
	static double grad[3];	// will be overwritten even called from different objects of this class
	grad[0] = ( q[idx] - q[im] ) / dx;	// 1 <= j <= Ny-1
	grad[1] = ( q[idx] - q[jm] ) / h[j];// 2 <= j <= Ny-1 (j==1 may not be valid)
	grad[2] = ( q[idx] - q[km] ) / dz;	// 1 <= j <= Ny-1
	return grad;
}


/****************************************/








// double* Scla::layerCC2UG(double *dst, double *src, int j1, int j0)
// /* interpolate j0 layer of src from cell center to U grid, result stored in j1 layer of dst */
// // note: pointer dst must NOT be the same as pointer src
// {
// 	int i, k;
// 	for (k=0; k<Nz; k++) {
// 	for (i=0; i<Nx; i++) {	dst[IDX(i,j1,k)] = 0.5 * ( src[IDX(i,j0,k)] + src[IDX(ima[i],j0,k)] );	}}
// 	return dst;
// }
// double* Scla::layerCC2VG(double *dst, double *src, int j1, int j0)
// /* interpolate j0,j0-1 layers of src from cell center to V grid, result stored in j1 layer of dst */
// // CAUTION: pointer dst must NOT be the same as pointer src
// {
// 	int i, k, jd = j0-1, ju = j0; // 1 <= j0 <= Ny
// 	for (k=0; k<Nz; k++) {
// 	for (i=0; i<Nx; i++) {
// 		dst[IDX(i,j1,k)] = ( src[IDX(i,jd,k)] * dy[ju] + src[IDX(i,ju,k)] * dy[jd] ) / (2.0*h[ju]);
// 	}}
// 	return dst;
// }
// double* Scla::layerCC2WG(double *dst, double *src, int j1, int j0)
// /* interpolate j0 layer of src from cell center to W grid, result stored in j1 layer of dst */
// // note: pointer dst must NOT be the same as pointer src
// {
// 	int i, k;
// 	for (k=0; k<Nz; k++) {
// 	for (i=0; i<Nx; i++) {	dst[IDX(i,j1,k)] = 0.5 * ( src[IDX(i,j0,k)] + src[IDX(i,j0,kma[k])] );	}}
// 	return dst;
// }


# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>


# include "Basic.h"

using namespace std;


/********** integration functions **********/

double Scla::yMeanU(int i, int k) const
{
	int j, ki = Nx * k + i;
	double integral = 0.0;
	for (j=1; j<Ny; j++) integral += this->q[Nxz * j + ki] * dy[j];
	return integral / Ly;
}
double Scla::yMeanV(int i, int k) const
{
	int j, ki = Nx * k + i;
	double integral = 0.0;
	for (j=1; j<=Ny; j++) integral += this->q[Nxz * j + ki] * (j==1 ? dy[1]/2 : j==Ny? dy[Ny-1]/2 : h[j]);
	return integral / Ly;
}
double Scla::layerMean(int j) const
{
	int k, i, kk, jj = Nxz * j;
	double integral = 0.0, dA = dx * dz, A = Lx * Lz;
	for (k=0; k<Nz; k++) { kk = jj + Nx * k;
	for (i=0; i<Nx; i++) { integral += this->q[kk + i] * dA;	}}
	return integral / A;
}
double Scla::bulkMeanU() const // only for U, W
{
	double integral = 0.0;
	for (int j=1; j<Ny; j++) integral += this->layerMean(j) * dy[j];
	return integral / Ly;
}
double Scla::bulkMeanV() const // only for V
{
	double integral = 0.0;
	for (int j=1; j<=Ny; j++) integral += this->layerMean(j) * (j==1 ? dy[1]/2 : j==Ny? dy[Ny-1]/2 : h[j]);
	return integral / Ly;
}


/********** interpolation functions **********/

void Scla::layerUG2CC(double *dst, int j)
/* interpolated quantity from U-grid to layer j of cell-center-grid */
{
	int i, k, kk;
	double *ql = this->lyrGet(j);
	for (k=0; k<Nz; k++) { kk = Nx * k;
	for (i=0; i<Nx; i++) { dst[kk+i] = 0.5 * (ql[kk+i] + ql[kk+ipa[i]]); }}
	
	if (j == 0 || j == Ny) { // interpolate virtual boundary to real boundary
		int ofst = (j==0 ? 1 : -1);
		double htmp = 0.5 / (j==0 ? h[1] : h[Ny]);
		ql = this->lyrGet(j + ofst);

		for (k=0; k<Nz; k++) { kk = Nx * k;
		for (i=0; i<Nx; i++) {
			dst[kk+i] = htmp * ( dst[kk+i] * dy[j+ofst] + .5 * (ql[kk+i] + ql[kk+ipa[i]]) * dy[j] );
		}}
	}
}

void Scla::layerWG2CC(double *dst, int j)
/* interpolated quantity from W-grid to layer j of cell-center-grid */
{
	int i, k, kk, kkf;
	double *ql = this->lyrGet(j);
	for (k=0; k<Nz; k++) { kk = Nx * k; kkf = Nx * kpa[k];
	for (i=0; i<Nx; i++) { dst[kk+i] = 0.5 * (ql[kk+i] + ql[kkf+i]); }}
	
	if (j == 0 || j == Ny) { // interpolate virtual boundary to real boundary
		int ofst = (j==0 ? 1 : -1);
		double htmp = 0.5 / (j==0 ? h[1] : h[Ny]);
		ql = this->lyrGet(j + ofst);

		for (k=0; k<Nz; k++) { kk = Nx * k; kkf = Nx * kpa[k];
		for (i=0; i<Nx; i++) {
			dst[kk+i] = htmp * ( dst[kk+i] * dy[j+ofst] + .5 * (ql[kk+i] + ql[kkf+i]) * dy[j] );
		}}
	}
}

void Scla::layerVG2CC(double *dst, int j)
/* interpolated quantity from V-grid to layer j of cell-center-grid */
{
	int i, k, kk;
	double *ql1 = this->lyrGet(j==0 ? 1 : j), *ql2 = this->lyrGet(j==Ny ? Ny : j+1);
	for (k=0; k<Nz; k++) { kk = Nx * k;
	for (i=0; i<Nx; i++) { dst[kk+i] = 0.5 * (ql1[kk+i] + ql2[kk+i]); }}
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


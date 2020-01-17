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
	double *q = blkGet(), integral = 0.0;
	for (j=1; j<Ny; j++)
		integral += q[Nxz * j + ki] * dy[j];
	return integral / Ly;
}
double Scla::yMeanV(int i, int k) const
{
	int j, ki = Nx * k + i;
	double *q = blkGet(), integral = 0.0;
	for (j=1; j<=Ny; j++)
		integral += q[Nxz * j + ki] * (j==1 ? dy[1]/2 : j==Ny ? dy[Ny-1]/2 : h[j]);
	return integral / Ly;
}
double Scla::layerMean(int j) const
{
	int k, i, kk, jj = Nxz * j;
	double *q = blkGet(), integral = 0.0, dA = dx * dz, A = Lx * Lz;
	for (k=0; k<Nz; k++) { kk = jj + Nx * k;
	for (i=0; i<Nx; i++) { integral += q[kk + i] * dA;	}}
	return integral / A;
}
double Scla::bulkMeanU() const // only for U, W
{
	double integral = 0.0;
	for (int j=1; j<Ny; j++)
		integral += layerMean(j) * dy[j];
	return integral / Ly;
}
double Scla::bulkMeanV() const // only for V
{
	double integral = 0.0;
	for (int j=1; j<=Ny; j++)
		integral += layerMean(j) * (j==1 ? dy[1]/2 : j==Ny? dy[Ny-1]/2 : h[j]);
	return integral / Ly;
}



/********** interpolation functions **********/

void Scla::layerUG2CC(double *dst, int j) const
/* interpolate quantity from U-grid to layer j of cell-center-grid */
{
	int i, k, kk; double *ql = lyrGet(j);
	for (k=0; k<Nz; k++) { kk = Nx * k;
	for (i=0; i<Nx; i++) { dst[kk+i] = 0.5 * (ql[kk+i] + ql[kk+ipa[i]]); }}
	
	if (j == 0 || j == Ny) { // interpolate virtual boundary to real boundary
		int ofst = (j==0 ? 1 : -1);
		double htmp = 0.5 / (j==0 ? h[1] : h[Ny]);
		ql = lyrGet(j + ofst);

		for (k=0; k<Nz; k++) { kk = Nx * k;
		for (i=0; i<Nx; i++) {
			dst[kk+i] = htmp * ( dst[kk+i] * dy[j+ofst] + .5 * (ql[kk+i] + ql[kk+ipa[i]]) * dy[j] );
		}}
	}
}

void Scla::layerWG2CC(double *dst, int j) const
/* interpolate quantity from W-grid to layer j of cell-center-grid */
{
	int i, k, kk, kkf; double *ql = lyrGet(j);
	for (k=0; k<Nz; k++) { kk = Nx * k; kkf = Nx * kpa[k];
	for (i=0; i<Nx; i++) { dst[kk+i] = 0.5 * (ql[kk+i] + ql[kkf+i]); }}
	
	if (j == 0 || j == Ny) { // interpolate virtual boundary to real boundary
		int ofst = (j==0 ? 1 : -1);
		double htmp = 0.5 / (j==0 ? h[1] : h[Ny]);
		ql = lyrGet(j + ofst);

		for (k=0; k<Nz; k++) { kk = Nx * k; kkf = Nx * kpa[k];
		for (i=0; i<Nx; i++) {
			dst[kk+i] = htmp * ( dst[kk+i] * dy[j+ofst] + .5 * (ql[kk+i] + ql[kkf+i]) * dy[j] );
		}}
	}
}

void Scla::layerVG2CC(double *dst, int j) const
/* interpolate quantity from V-grid to layer j of cell-center-grid */
{
	int i, k, kk;
	double *ql1 = lyrGet(j==0 ? 1 : j), *ql2 = lyrGet(j==Ny ? Ny : j+1);
	for (k=0; k<Nz; k++) { kk = Nx * k;
	for (i=0; i<Nx; i++) { dst[kk+i] = 0.5 * (ql1[kk+i] + ql2[kk+i]); }}
}



void Scla::CC2EG(double *dst1, double *dst2, double *dst3) const
{
	int i, j, k, idx, im, jm, km, imjm, jmkm, imkm; double *q = blkGet();
	for (j=0; j<=Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k); im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
		imjm = IDX(ima[i],j-1,k); jmkm = IDX(i,j-1,kma[k]); imkm = IDX(ima[i],j,kma[k]);

		dst2[idx] = .25 * ( q[idx] + q[im] + q[km] + q[imkm] );                        if (j>0)
		dst1[idx] = .25/h[j] * ( (q[idx]+q[km]) * dy[j-1] + (q[jm]+q[jmkm]) * dy[j] ); if (j>0)
		dst3[idx] = .25/h[j] * ( (q[idx]+q[im]) * dy[j-1] + (q[jm]+q[imjm]) * dy[j] );
	}}}
}


void Scla::layerCC2XE(double *dst, int j) const // j = 1 ~ Ny
/* interpolate quantity from cell-center to layer j of X edge */
{
	int i, k, idx, km; double *ql1 = lyrGet(j-1), *ql2 = lyrGet(j);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = Nx * k + i;
		km = Nx * kma[k] + i;
		dst[idx] = .25/h[j] * ( (ql2[idx]+ql2[km]) * dy[j-1] + (ql1[idx]+ql1[km]) * dy[j] );
	}}
}

void Scla::layerCC2ZE(double *dst, int j) const // j = 1 ~ Ny
/* interpolate quantity from cell-center to layer j of Z edge */
{
	int i, k, idx, im; double *ql1 = lyrGet(j-1), *ql2 = lyrGet(j);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = Nx * k + i;
		im = Nx * k + ima[i];
		dst[idx] = .25/h[j] * ( (ql2[idx]+ql2[im]) * dy[j-1] + (ql1[idx]+ql1[im]) * dy[j] );
	}}
}

void Scla::layerCC2YE(double *dst, int j) const // j = 0 ~ Ny
/* interpolate quantity from cell-center to layer j of Y edge */
{
	int i, k, idx, im, km, imkm; double *ql = lyrGet(j);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = Nx * k + i;
		im = Nx * k + ima[i];
		km = Nx * kma[k] + i;
		imkm = Nx * kma[k] + ima[i];
		dst[idx] = .25 * ( ql[idx] + ql[im] + ql[km] + ql[imkm] );
	}}
}



/********** differentiation operators **********/

double* Scla::gradient(int i, int j, int k) const
/* compute gradient of a scalar field at cell center to grid points corresponding to U,V,W */
// CAUTION: avoid successive calling to this function, because the static return variable will be overwritten every time
{
	int idx= IDX(i,j,k);
	int im = IDX(ima[i],j,k);
	int jm = IDX(i,j-1,k);
	int km = IDX(i,j,kma[k]);
	double *q = blkGet();
	static double grad[3];	// will be overwritten even called from different objects of this class
	grad[0] = ( q[idx] - q[im] ) / dx;	// 1 <= j <= Ny-1
	grad[1] = ( q[idx] - q[jm] ) / h[j];// 2 <= j <= Ny-1 (j==1 may not be valid)
	grad[2] = ( q[idx] - q[km] ) / dz;	// 1 <= j <= Ny-1
	return grad;
}


/****************************************/










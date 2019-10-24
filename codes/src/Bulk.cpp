# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>


# include "Basic.h"
typedef fftw_complex fcmplx;

using namespace std;


Bulk::Bulk(int n1, int n2, int n3, bool inift):
nx(n1),
ny(n2),
nz(n3),
nxz(n1*n3),
nxc(nx/2+1),
nxr(2*nxc),
nxzc(nz*nxc),
nxzr(nz*nxr)
{
	q = new double [nxz*ny];

	if (inift) {
		fq = new double [nxzr*ny];
		frcs = new fftw_plan [ny];
		fcrs = new fftw_plan [ny];
		for (int j=0; j<ny; j++) {	// note: the effective domain of DP actually starts from the second layer, since the wall layers do not participate in computation
			double *r = (double*) &(q[nxz * j]);
			fcmplx *c = (fcmplx*) &(fq[nxzr * j]);
			frcs[j] = fftw_plan_dft_r2c_2d(nz, nx, r, c, FFTW_MEASURE);
			fcrs[j] = fftw_plan_dft_c2r_2d(nz, nx, c, r, FFTW_MEASURE); // note: the inverse transform (complex to real) has the side-effect of overwriting its input array
		}	
	}
	else {
		fq = NULL;
		frcs = NULL;
		fcrs = NULL;
	}
}

void Bulk::freeall()
{
	delete [] q;
	if (fq) {
		delete [] fq;
		for (int j=0; j<ny; j++) {
			fftw_destroy_plan(frcs[j]);
			fftw_destroy_plan(fcrs[j]);
		}
		delete [] frcs;
		delete [] fcrs;
	}
}



double* Bulk::fft()
{
	for (int j=0; j<ny; j++) fftw_execute(frcs[j]);
	return fq;
}

double* Bulk::ifft()
{
	for (int j=0; j<ny; j++) fftw_execute(fcrs[j]);
	this->bulkMlt(1.0/nxz);
	return q;
}


/***** convinent operations for whole arrays *****/
double* Bulk::layerGet(int j)
{ return (double*) &(q[nxz * j]); }

double* Bulk::layerCpy(double *src, int j1, int j0)
{
	memcpy(& q[nxz * j1], & src[nxz * j0], sizeof(double) * nxz);
	return q;
}
double* Bulk::layerSet(double a, int j)
{
	int k, i, kk, jj = nxz * j;
	for (k=0; k<nz; k++) { kk = jj + nx * k;
	for (i=0; i<nx; i++) { q[kk + i] = a;	}}
	return q;
}
double* Bulk::layerAdd(double a, int j)
{
	int k, i, kk, jj = nxz * j;
	for (k=0; k<nz; k++) { kk = jj + nx * k;
	for (i=0; i<nx; i++) { q[kk + i] += a;	}}
	return q;
}
double* Bulk::layerMlt(double a, int j)
{
	int k, i, kk, jj = nxz * j;
	for (k=0; k<nz; k++) { kk = jj + nx * k;
	for (i=0; i<nx; i++) { q[kk + i] *= a;	}}
	return q;
}
double* Bulk::layersAdd(double *src, int j1, int j0)
{	// CAUTION: if src is self, different layers may interfere
	int k, i, kk0, kk1, jj0 = nxz * j0, jj1 = nxz * j1;
	for (k=0; k<nz; k++) { kk0 = jj0 + nx * k; kk1 = jj1 + nx * k;
	for (i=0; i<nx; i++) { q[kk1 + i] += src[kk0 + i];	}}
	return q;
}
double* Bulk::layersMlt(double *src, int j1, int j0)
{	// CAUTION: if src is self, different layers may interfere
	int k, i, kk0, kk1, jj0 = nxz * j0, jj1 = nxz * j1;
	for (k=0; k<nz; k++) { kk0 = jj0 + nx * k; kk1 = jj1 + nx * k;
	for (i=0; i<nx; i++) { q[kk1 + i] *= src[kk0 + i];	}}
	return q;
}


double* Bulk::bulkCpy(double *src)
{	// note: bulk functions will change the boundary
	memcpy(q, src, sizeof(double) * nxz * ny);
	return q;
}
double* Bulk::bulkSet(double a)
{	// note: bulk functions will change the boundary
	for (int j=0; j<ny; j++)	this->layerSet(a, j);
	return q;
}
double* Bulk::bulkAdd(double a)
{	// note: bulk functions will change the boundary
	for (int j=0; j<ny; j++)	this->layerAdd(a, j);
	return q;
}
double* Bulk::bulkMlt(double a)
{	// note: bulk functions will change the boundary
	for (int j=0; j<ny; j++)	this->layerMlt(a, j);
	return q;
}


/***** file IO *****/

void Bulk::fileIO(char *path, char *name, char mode)
/* read & write field from & to binary files */
{
	FILE *fp;
	char str[1024];

	sprintf(str, "%s%s.bin", path, name);

	fp = fopen( str, (mode == 'w' ? "wb" : "rb") );

		// write domain information at the beginning
		if (mode == 'w') {
			fwrite(& nx, sizeof(int), 1, fp);
			fwrite(& ny, sizeof(int), 1, fp);
			fwrite(& nz, sizeof(int), 1, fp);
		}
		// data begin after the info section
		fseek(fp, sizeof(double) * nxz, SEEK_SET);
		if (mode == 'w') fwrite(q, sizeof(double) * nxz, ny, fp);
		if (mode == 'r') fread (q, sizeof(double) * nxz, ny, fp);

	fclose(fp);
}

void Bulk::debug_AsciiOutput(char *path, char *name, int j1, int j2)
/* write the fields in ascii files for check */
{
	FILE *fp;
	char str[1024];
	int i, j, k;

	sprintf(str, "%s%s.txt", path, name);
	fp = fopen(str, "w");
	for (j=j1;j<j2; j++) { fputc('\n', fp);
	for (k=0; k<nz; k++) { fputc('\n', fp);
	for (i=0; i<nx; i++) {
		sprintf(str, "%.6f\t", this->id(i,j,k));
		fputs(str, fp);
	}}}
	fclose(fp);

	if (this->fq) {

		sprintf(str, "%s%s_FR.txt", path, name);
		fp = fopen(str, "w");
		for (j=j1;j<j2; j++) { fputc('\n', fp);
		for (k=0; k<nz; k++) { fputc('\n', fp);
		for (i=0; i<nxc;i++) {
			sprintf(str, "%.6f\t", this->idf(2*i,j,k));
			fputs(str, fp);
		}}}
		fclose(fp);

		sprintf(str, "%s%s_FI.txt", path, name);
		fp = fopen(str, "w");
		for (j=j1;j<j2; j++) { fputc('\n', fp);
		for (k=0; k<nz; k++) { fputc('\n', fp);
		for (i=0; i<nxc;i++) {
			sprintf(str, "%.6f\t", this->idf(2*i+1,j,k));
			fputs(str, fp);
		}}}
		fclose(fp);

	}
}





# include <iostream>
# include <stdlib.h>

# include "Matrix.h"

using namespace std;



Matrix::Matrix(int rank):	n(rank)
{
	// l = (double*) malloc(sizeof(double) * n);	// array l is not needed since it can be replaced by a single number l in execution
	u = (double*) malloc(sizeof(double) * n);
	rho = (double*) malloc(sizeof(double) * n);
	sig = (double*) malloc(sizeof(double) * n);
}

Matrix::~Matrix()
{
	// free(l);
	free(u);
	free(rho);
	free(sig);
}

void Matrix::tdma(double *a, double *b, double *c, double *d)
/* Tridiagonal matrix algorithm, result returned in d. (a:1~n-1, b:0~n-1, c:0~n-2, d:0~n-1) */
/* example n = 5:
|b0 c0 0  0  0 |   |d0|
|a1 b1 c1 0  0 |   |d1|
|0  a2 b2 c2 0 | = |d2|
|0  0  a3 b3 c3|   |d3|
|0  0  0  a4 b4|   |d4|	*/
{
	int i;
	double l;

	u[0] = b[0];
	for (i=1; i<n; i++) {
		l = a[i] / u[i-1];
		u[i] = b[i] - l * c[i-1];
		d[i] -= l * d[i-1];	// solve intermediate variable y and store in d
	}

	d[n-1] /= u[n-1];
	for (i=n-2; i>=0; i--) {
		d[i] = ( d[i] - c[i] * d[i+1] ) / u[i];
	}
}

void Matrix::ctdma(double *a, double *b, double *c, double *d)
/* Periodic tridiagonal matrix algorithm, result returned in d. (a:0~n-1, b:0~n-1, c:0~n-1, d:0~n-1) */
/* example n = 5:
|b0 c0 0  0  a0|   |d0|
|a1 b1 c1 0  0 |   |d1|
|0  a2 b2 c2 0 | = |d2|
|0  0  a3 b3 c3|   |d3|
|c4 0  0  a4 b4|   |d4|	*/
{
	int i;
	double l, sum1 = 0.0, sum2 = 0.0;

	u[0] = b[0];
	rho[0] = a[0];
	sig[0] = c[n-1] / u[0];

	for (i=1; i<n; i++) {
		l = a[i] / u[i-1];
		u[i] = b[i] - l * c[i-1];
		rho[i] = - l * rho[i-1]; // rho[n-1] is useless
		sig[i] = - c[i-1] * sig[i-1] / u[i]; // sig[n-1] is useless
		d[i] -= l * d[i-1];	// solve intermediate variable y and store in d

		sum1 += rho[i-1] * sig[i-1];
		sum2 += sig[i-1] * d[i-1];
	}

	i = n-1;
	u[i] -= l*rho[i-1] + c[i-1]*sig[i-1] + sum1;
	d[i] -= sum2;
	d[i] /= u[i];

	for (i=n-2; i>=0; i--) {
		d[i] = ( d[i] - c[i] * d[i+1] - rho[i] * d[n-1] ) / u[i];
	}
}


// this function does not belong to the Matrix class
double linearInterp(int n, double *x0, double *y0, double x)
{
	int i;
	for (i=0; i<n; i++) { if (x0[i] > x) break; }
	return	(i==0) ? y0[0]
		:	(i==n) ? y0[n-1]
		:	(	( x0[i] - x ) * y0[i-1]
			+	( x - x0[i-1] ) * y0[i]	) / ( x0[i] - x0[i-1] );
}



// # define DEBUG	// g++ -lfftw3 -lm -I include src/Matrix.cpp

# ifdef DEBUG

void tdma_validation()
{
	// tdma
	double a[4] = {1,2,3,4};
	double b[4] = {5,2,3,4};
	double c[4] = {1,2,3,4};
	double d[4] = {1,2,3,4};

	cout << "\ntdma validation:" << endl;
	tdma(4,a,b,c,d);
	printf("%.4f\n%.4f\n%.4f\n%.4f\n", d[0], d[1], d[2], d[3]);

	// ctdma
	double e[5] = {5,4,3,2,1};
	double f[5] = {5,2,3,4,5};
	double g[5] = {1,2,3,4,6};
	double h[5] = {5,4,3,2,1};

	cout << "\nctdma validation:" << endl;
	ctdma(5,e,f,g,h);
	printf("%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n", h[0], h[1], h[2], h[3], h[4]);
}

int main()
{
	tdma_validation();
}

# endif




# ifdef MATRIX_NO_CLASS

void tdma(int n, double *a, double *b, double *c, double *d)
/* Tridiagonal matrix algorithm, result returned in d. (a:1~n-1, b:0~n-1, c:0~n-2, d:0~n-1) */
/* example n = 5:
|b0 c0 0  0  0 |   |d0|
|a1 b1 c1 0  0 |   |d1|
|0  a2 b2 c2 0 | = |d2|
|0  0  a3 b3 c3|   |d3|
|0  0  0  a4 b4|   |d4|	*/
{
	int i;
	double l;
	double *u = (double*) malloc(sizeof(double) * n);

	u[0] = b[0];
	for (i=1; i<n; i++) {
		l = a[i] / u[i-1];
		u[i] = b[i] - l * c[i-1];
		d[i] -= l * d[i-1];	// solve intermediate variable y and store in d
	}

	d[n-1] /= u[n-1];
	for (i=n-2; i>=0; i--) {
		d[i] = ( d[i] - c[i] * d[i+1] ) / u[i];
	}

	free(u);
}

void ctdma(int n, double *a, double *b, double *c, double *d)
/* Periodic tridiagonal matrix algorithm, result returned in d. (a:0~n-1, b:0~n-1, c:0~n-1, d:0~n-1) */
/* example n = 5:
|b0 c0 0  0  a0|   |d0|
|a1 b1 c1 0  0 |   |d1|
|0  a2 b2 c2 0 | = |d2|
|0  0  a3 b3 c3|   |d3|
|c4 0  0  a4 b4|   |d4|	*/
{
	int i;
	double l, sum1 = 0.0, sum2 = 0.0;
	double *u = (double*) malloc(sizeof(double) * n);
	double *rho = (double*) malloc(sizeof(double) * n);
	double *sig = (double*) malloc(sizeof(double) * n);

	u[0] = b[0];
	rho[0] = a[0];
	sig[0] = c[n-1] / u[0];

	for (i=1; i<n; i++) {
		l = a[i] / u[i-1];
		u[i] = b[i] - l * c[i-1];
		rho[i] = - l * rho[i-1]; // rho[n-1] is useless
		sig[i] = - c[i-1] * sig[i-1] / u[i]; // sig[n-1] is useless
		d[i] -= l * d[i-1];	// solve intermediate variable y and store in d

		sum1 += rho[i-1] * sig[i-1];
		sum2 += sig[i-1] * d[i-1];
	}

	i = n-1;
	u[i] -= l*rho[i-1] + c[i-1]*sig[i-1] + sum1;
	d[i] -= sum2;
	d[i] /= u[i];

	for (i=n-2; i>=0; i--) {
		d[i] = ( d[i] - c[i] * d[i+1] - rho[i] * d[n-1] ) / u[i];
	}

	free(u);
	free(rho);
	free(sig);
}

# endif








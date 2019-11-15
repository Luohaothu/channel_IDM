# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>


# include "Basic.h"

using namespace std;




double Vctr::divergence(int i, int j, int k)
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

double Vctr::convection(int i, int j, int k)
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

double* Vctr::strainrate(int i, int j, int k)
/* compute the strain rate tensor of a vector field at the center of cell (i,j,k) */
// CAUTION: avoid successive calling to this function, because the static return variable will be overwritten every time
{
	int idx = IDX(i,j,k);	// 1 <= j <= Ny-1
	int ip = IDX(ipa[i],j,k), jp = IDX(i,j+1,k), kp = IDX(i,j,kpa[k]);
	int im = IDX(ima[i],j,k), jm = IDX(i,j-1,k), km = IDX(i,j,kma[k]);
	int ipjp = IDX(ipa[i],j+1,k), jpkp = IDX(i,j+1,kpa[k]), ipkp = IDX(ipa[i],j,kpa[k]);
	int ipjm = IDX(ipa[i],j-1,k), jpkm = IDX(i,j+1,kma[k]), imkp = IDX(ima[i],j,kpa[k]);
	int imjp = IDX(ima[i],j+1,k), jmkp = IDX(i,j-1,kpa[k]), ipkm = IDX(ipa[i],j,kma[k]);
	double u1, u2, v1, v2, w1, w2;
	static double sr[6];	// will be overwritten even called from different objects of this class

	// S_ii
	sr[0] = ( u[ip] - u[idx] ) / dx;
	sr[1] = ( v[jp] - v[idx] ) / dy[j];
	sr[2] = ( w[kp] - w[idx] ) / dz;
	// S_12
	v2 = 0.25 * ( v[idx] + v[ip] + v[jp] + v[ipjp] );
	v1 = 0.25 * ( v[idx] + v[im] + v[jp] + v[imjp] );
	u2 = ( (u[idx]+u[ip]) * dy[j+1] + (u[jp]+u[ipjp]) * dy[j] ) / (4.0*h[j+1]);
	u1 = ( (u[idx]+u[ip]) * dy[j-1] + (u[jm]+u[ipjm]) * dy[j] ) / (4.0*h[j]);
	sr[3] = 0.5 * ( (v2-v1) / dx + (u2-u1) / dy[j] );
	// S_23
	w2 = ( (w[idx]+w[kp]) * dy[j+1] + (w[jp]+w[jpkp]) * dy[j] ) / (4.0*h[j+1]);
	w1 = ( (w[idx]+w[kp]) * dy[j-1] + (w[jm]+w[jmkp]) * dy[j] ) / (4.0*h[j]);
	v2 = 0.25 * ( v[idx] + v[jp] + v[kp] + v[jpkp] );
	v1 = 0.25 * ( v[idx] + v[jp] + v[km] + v[jpkm] );
	sr[4] = 0.5 * ( (w2-w1) / dy[j] + (v2-v1) / dz );
	// S_13
	w2 = 0.25 * ( w[idx] + w[ip] + w[kp] + w[ipkp] );
	w1 = 0.25 * ( w[idx] + w[im] + w[kp] + w[imkp] );
	u2 = 0.25 * ( u[idx] + u[ip] + u[kp] + u[ipkp] );
	u1 = 0.25 * ( u[idx] + u[ip] + u[km] + u[ipkm] );
	sr[5] = 0.5 * ( (w2-w1) / dx + (u2-u1) / dz );

	return sr;
}

double* Vctr::gradient(int i, int j, int k)
/* compute the gradient tensor of a vector field at the center of cell (i,j,k) */
// CAUTION: avoid successive calling to this function, because the static return variable will be overwritten every time
{
	int idx = IDX(i,j,k);
	int ip = IDX(ipa[i],j,k), jp = IDX(i,j+1,k), kp = IDX(i,j,kpa[k]);
	int im = IDX(ima[i],j,k), jm = IDX(i,j-1,k), km = IDX(i,j,kma[k]);
	int ipjp = IDX(ipa[i],j+1,k), jpkp = IDX(i,j+1,kpa[k]), ipkp = IDX(ipa[i],j,kpa[k]);
	int ipjm = IDX(ipa[i],j-1,k), jpkm = IDX(i,j+1,kma[k]), imkp = IDX(ima[i],j,kpa[k]);
	int imjp = IDX(ima[i],j+1,k), jmkp = IDX(i,j-1,kpa[k]), ipkm = IDX(ipa[i],j,kma[k]);
	double u1, u2, v1, v2, w1, w2;
	static double gr[9];	// will be overwritten even called from different objects of this class

	// a11
	gr[0] = ( u[ip] - u[idx] ) / dx;
	// a12
	v2 = 0.25 * ( v[idx] + v[ip] + v[jp] + v[ipjp] );
	v1 = 0.25 * ( v[idx] + v[im] + v[jp] + v[imjp] );
	gr[1] = (v2-v1) / dx;
	// a13
	w2 = 0.25 * ( w[idx] + w[ip] + w[kp] + w[ipkp] );
	w1 = 0.25 * ( w[idx] + w[im] + w[kp] + w[imkp] );
	gr[2] = (w2-w1) / dx;

	// a21
	u2 = ( (u[idx]+u[ip]) * dy[j+1] + (u[jp]+u[ipjp]) * dy[j] ) / (4.0*h[j+1]);
	u1 = ( (u[idx]+u[ip]) * dy[j-1] + (u[jm]+u[ipjm]) * dy[j] ) / (4.0*h[j]);
	gr[3] = (u2-u1) / dy[j];
	// a22
	gr[4] = ( v[jp] - v[idx] ) / dy[j];
	// a23
	w2 = ( (w[idx]+w[kp]) * dy[j+1] + (w[jp]+w[jpkp]) * dy[j] ) / (4.0*h[j+1]);
	w1 = ( (w[idx]+w[kp]) * dy[j-1] + (w[jm]+w[jmkp]) * dy[j] ) / (4.0*h[j]);
	gr[5] = (w2-w1) / dy[j];
	
	// a31
	u2 = 0.25 * ( u[idx] + u[ip] + u[kp] + u[ipkp] );
	u1 = 0.25 * ( u[idx] + u[ip] + u[km] + u[ipkm] );
	gr[6] = (u2-u1) / dz;
	// a32
	v2 = 0.25 * ( v[idx] + v[jp] + v[kp] + v[jpkp] );
	v1 = 0.25 * ( v[idx] + v[jp] + v[km] + v[jpkm] );
	gr[7] = (v2-v1) / dz;
	// a33
	gr[8] = ( w[kp] - w[idx] ) / dz;

	return gr;
}







// double Vctr::divergence(int i, int j, int k)
// /* compute divergence of a vector field at the center of cell (i,j,k) */
// {
// 	return	(u.id(ipa[i],j,k) - u.id(i,j,k)) / u.dx
// 		+	(v.id(i,j+1,k)    - v.id(i,j,k)) / v.dy[j]	// 1 <= j <= Ny-1
// 		+	(w.id(i,j,kpa[k]) - w.id(i,j,k)) / w.dz;
// }

// double Vctr::convection(int i, int j, int k)
// /* compute convection coefficient of a vector field at the center of cell (i,j,k) */
// {
// 	return	0.5 * (u.id(ipa[i],j,k) + u.id(i,j,k)) / u.dx
// 		+	0.5 * (v.id(i,j+1,k)    + v.id(i,j,k)) / v.dy[j]	// 1 <= j <= Ny-1
// 		+	0.5 * (w.id(i,j,kpa[k]) + w.id(i,j,k)) / w.dz;
// }
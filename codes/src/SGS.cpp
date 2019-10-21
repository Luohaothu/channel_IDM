# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>


# include "SGS.h"


void SGS::viscalc(double *nu, double *U[3], double Re, class Mesh *pmesh)
{
	int i, j, k;
	this->smargorinsky(nu, U[0], U[1], U[2], Re, pmesh);
	for (j=0; j<=Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {	nu[IDX(i,j,k)] += 1.0/Re;	}}}
}

void SGS::smargorinsky(double *ev, double *u, double *v, double *w, double Re, class Mesh *pmesh)
{
	int i, j, k, idx;
	double *sr, sra, dlt, dmp;

	// calculate eddy viscosity at cell centers
	for (j=1; j<Ny; j++) {

		dlt = pow( pmesh->dx * pmesh->dy[j] * pmesh->dz, 1.0/3.0 ); // filter size
		dmp = 1.0 - exp( -pmesh->y[j] * Re / 22.0 / 25.0 ); // Van Driest damping

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			sr = pmesh->strainrate(u, v, w, i, j, k);
			sra = sqrt( 2.0 * (
					sr[0]*sr[0] + sr[1]*sr[1] + sr[2]*sr[2]
				+ 2.0 * ( sr[3]*sr[3] + sr[4]*sr[4] + sr[5]*sr[5] )
			) );

			ev[IDX(i,j,k)] = cs * pow(dlt * dmp, 2.0) * sra;
		}}
	}
	// equate the boundary eddy viscosity to the first layer off wall
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		ev[IDX(i,0,k)] = ev[IDX(i,1,k)];
		ev[IDX(i,Ny,k)] = ev[IDX(i,Ny-1,k)];
	}}
}










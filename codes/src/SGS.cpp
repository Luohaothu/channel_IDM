# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>


# include "SGS.h"


SGS::SGS(int *dim)
Nx(dim[0]), Ny(dim[1]), Nz(dim[2]), Nxz(dim[0]*dim[2])
{
	ev = new double [Nxz * (Ny+1)];
	evi = new double [Nxz * (Ny+1)];
	evj = new double [Nxz * (Ny+1)];
	evk = new double [Nxz * (Ny+1)];
}


void SGS::smargorinsky(double *u, double *v, double *w, double Re, class Mesh *pmesh)
{
	int i, j, k, idx;
	double *sr, sra, dlt, dmp, cs = 0.01;

	// calculate eddy viscosity at cell centers
	for (j=1; j<Ny; j++) {

		dlt = pow( pmesh->dx * pmesh->dy[j] * pmesh->dz, 1.0/3.0 ); // filter size
		dmp = 1.0 - exp( -pmesh->y[j] * Re / 22.0 / 25.0 ); // Van Driest damping

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			idx = IDX(i,j,k);

			sr = pmesh->stainrate(u, v, w, i, j, k);
			sra = sqrt( 2.0 * (
					sr[0]*sr[0] + sr[1]*sr[1] + sr[2]*sr[2]
				+ 2.0 * ( sr[3]*sr[3] + sr[4]*sr[4] + sr[5]*sr[5] )
			) );

			ev[idx] = cs * pow(dlt * dmp, 2.0) * sra;
		}}
	}
	// equate the boundary eddy viscosity to the nearest grid
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		ev[IDX(i,0,k)] = ev[IDX(i,1,k)];
		ev[IDX(i,Ny,k)] = ev[IDX(i,Ny-1,k)];
	}}
	
	// interpolate to axes
	for (j=0; j<=Ny; j++) {
	for (i=0; i<Nx; i++) {
	for (k=0; k<Nz; k++) {

	}}}
}
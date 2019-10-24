# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>


# include "SGS.h"


void SGS::evcalc(Scla &EV, Vctr &U, double Re, double Cs)
{
	this->Cs = Cs;
	this->smargorinsky(EV, U, Re);
}



void SGS::smargorinsky(Scla &EV, Vctr &U, double Re)
{
	int i, j, k, idx;
	double *sr, sra, dlt, dmp;

	// calculate eddy viscosity at cell centers
	for (j=1; j<Ny; j++) {

		dlt = pow( dvol[j], 1.0/3.0 ); // filter size
		dmp = 1.0 - exp( -y[j] * Re / 22.0 / 25.0 ); // Van Driest damping

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			sr = U.strainrate(i, j, k);
			sra = sqrt( 2.0 * (
					sr[0]*sr[0] + sr[1]*sr[1] + sr[2]*sr[2]
				+ 2.0 * ( sr[3]*sr[3] + sr[4]*sr[4] + sr[5]*sr[5] )
			) );

			EV.id(i,j,k) = Cs * pow(dlt * dmp, 2.0) * sra;
		}}
	}
	// equate the boundary eddy viscosity to the first layer off wall
	EV.layerCpy(EV, 0 , 1);
	EV.layerCpy(EV, Ny, Ny-1);

}










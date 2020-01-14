# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <cmath>

# include "Solver.h"
# include "Interp.h"

using namespace std;



void Solver::getbc()
/* set boundary conditions UBC (not yet applied to velocity field) */
{
	field.UBC[1] = 0.;
	field.UBC[2] = 0.;
	field.UBC[3] = 0.;
}

void Solver::getbc(Vctr &U0)
{
	Interp(U0[1], field.UBC[1]).bulkFilter('U');
	Interp(U0[2], field.UBC[2]).bulkFilter('X');
	Interp(U0[3], field.UBC[3]).bulkFilter('U');
}

void Solver::getnu(double Re, int bftype)
{
	switch (bftype) {
		case 2: sgs.smargorinsky(0.18, Re); field.NU += 1./Re; break;
		case 3: sgs.dynamicsmarg();         field.NU += 1./Re; break;
		case 4: sgs.dynamicvreman(Re);      field.NU += 1./Re; break;
		default:                            field.NU  = 1./Re;
	}
}

void Solver::getfb()
{
	field.FB[1] = 0.;
	field.FB[2] = 0.;
	field.FB[3] = 0.;
}

void Solver::removeSpanMean()
{
	int i, j, k, n;
	double qm, *qsm = new double [mesh.Nx];
	Bulk *bks[4] = {&field.U[1], &field.U[2], &field.UH[1], &field.UH[2]};

	for (n=0; n<4; n++) {
	for (j=1; j<mesh.Ny; j++) {
		if (j==1 && (n==1 || n==3)) continue; // j==1 layer for V & VH is boundary, skip

		for (i=0; i<mesh.Nx; i++) { qsm[i] = 0.0;
		for (k=0; k<mesh.Nz; k++) { qsm[i] += bks[n]->id(i,j,k) / mesh.Nz; }}

		qm = 0.0;
		for (i=0; i<mesh.Nx; i++)   qm += qsm[i] / mesh.Nx;

		for (k=0; k<mesh.Nz; k++) {
		for (i=0; i<mesh.Nx; i++) { bks[n]->id(i,j,k) -= qsm[i] - qm; }}
	}}

	delete [] qsm;
}
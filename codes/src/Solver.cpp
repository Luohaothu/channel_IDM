# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <cmath>

# include "Solver.h"
# include "Interp.h"

using namespace std;




void Solver::removeSpanMean()
{
	int i, j, k, n;
	double qm, *qsm = new double [mesh.Nx];
	Bulk *bks[4] = {&FLD.V[1], &FLD.V[2], &FLDH.V[1], &FLDH.V[2]};

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


void Solver::debug_Output(int tstep)
{
	char path[1024] = "debug/";
	{
		char names[4][32] = {"U", "V", "W", "P"};
		FLD.V[1].debug_AsciiOutput(path, names[0], 0, mesh.Ny+1);
		FLD.V[2].debug_AsciiOutput(path, names[1], 1, mesh.Ny+1);
		FLD.V[3].debug_AsciiOutput(path, names[2], 0, mesh.Ny+1);
		   FLD.S.debug_AsciiOutput(path, names[3], 0, mesh.Ny+1);
	}{
		char names[4][32] = {"UH", "VH", "WH", "DP"};
		FLDH.V[1].debug_AsciiOutput(path, names[0], 0, mesh.Ny+1);
		FLDH.V[2].debug_AsciiOutput(path, names[1], 1, mesh.Ny+1);
		FLDH.V[3].debug_AsciiOutput(path, names[2], 0, mesh.Ny+1);
		   FLDH.S.debug_AsciiOutput(path, names[3], 0, mesh.Ny+1);
	}{
		char names[4][32] = {"NUX", "NUY", "NUZ", "NU"};
		VIS.V[1].debug_AsciiOutput(path, names[0], 0, mesh.Ny+1);
		VIS.V[2].debug_AsciiOutput(path, names[1], 0, mesh.Ny+1);
		VIS.V[3].debug_AsciiOutput(path, names[2], 0, mesh.Ny+1);
		   VIS.S.debug_AsciiOutput(path, names[3], 0, mesh.Ny+1);
	}{
		char names[4][32] = {"UBC", "VBC", "WBC", "PBC"};
		BC.V[1].debug_AsciiOutput(path, names[0], 0, 2);
		BC.V[2].debug_AsciiOutput(path, names[1], 0, 2);
		BC.V[3].debug_AsciiOutput(path, names[2], 0, 2);
		   BC.S.debug_AsciiOutput(path, names[3], 0, 2);
	}{
		char names[3][32] = {"FBX", "FBY", "FBZ"};
		FB[1].debug_AsciiOutput(path, names[0], 0, mesh.Ny+1);
		FB[2].debug_AsciiOutput(path, names[1], 0, mesh.Ny+1);
		FB[3].debug_AsciiOutput(path, names[2], 0, mesh.Ny+1);
	}
}






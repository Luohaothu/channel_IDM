# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>

# include "Interp.h"
# include "SGS.h"


// void SGS::evcalc(Scla &EV, Vctr &U, double Re, double Cs)
// {
// 	this->Cs = Cs;

// 	if (Re && Cs) this->smargorinsky(EV, U, Re);
// 	else if (!(Re || Cs)) this->dynamicsmarg(EV, U);
// }



void SGS::smargorinsky(Scla &EV, Vctr &U, double Cs, double Re)
{
	int i, j, k;
	double *sr, sra, dlt, dmp;

	// calculate delta_nu
	Scla &U1 = U.com1;
	double tauw = (0.5/Re) * (
		(U1.layerMean(1) - U1.layerMean(0))    / h[1]
	-	(U1.layerMean(Ny)- U1.layerMean(Ny-1)) / h[Ny] );
	double Ret = Re * sqrt(tauw);	// Ly is taken for 2.0 as default

	// calculate eddy viscosity at cell centers
	for (j=1; j<Ny; j++) {

		dlt = pow( dvol[j], 1.0/3.0 ); // filter size
		dmp = 1.0 - exp( (fabs(1.-yc[j]) - 1.) * Ret / 25.0 ); // Van Driest damping: 1-exp(-y^+/A^+), with A^+ = 25

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			sr = U.strainrate(i,j,k);
			sra = sqrt( 2.0 * (
					sr[0]*sr[0] + sr[1]*sr[1] + sr[2]*sr[2]
				+ 2.0 * ( sr[3]*sr[3] + sr[4]*sr[4] + sr[5]*sr[5] )
			) );

			EV.id(i,j,k) = pow(Cs * dlt * dmp, 2.0) * sra;
		}}
	}
	// equate the boundary eddy viscosity to the first layer off wall
	EV.layerCpy(EV, 0 , 1);
	EV.layerCpy(EV, Ny, Ny-1);

}


void SGS::dynamicsmarg(Scla &EV, Vctr &U)
{
	int i, j, k;
	double *sr, sra, dlt1, dlt2, u, v, w;

	Mesh mesh0(Nx,0,Nz,Lx,0,Lz);
	Scla S11(mesh0), S22(mesh0), S33(mesh0), S12(mesh0), S23(mesh0), S13(mesh0); 
	Scla M11(mesh0), M22(mesh0), M33(mesh0), M12(mesh0), M23(mesh0), M13(mesh0);
	Scla L11(mesh0), L22(mesh0), L33(mesh0), L12(mesh0), L23(mesh0), L13(mesh0); 
	Scla &UC1 = L11, &UC2 = L22, &UC3 = L33, &UF1 = S12, &UF2 = S23, &UF3 = S13;

	Interp flt1 (UC1,UF1), flt2 (UC2,UF2), flt3 (UC3,UF3);
	Interp fsm11(S11,M11), fsm22(S22,M22), fsm33(S33,M33), \
		   fsm12(S12,M12), fsm23(S23,M23), fsm13(S13,M13);
	Interp fsl11(S11,L11), fsl22(S22,L22), fsl33(S33,L33), \
		   fsl12(S12,L12), fsl23(S23,L23), fsl13(S13,L13);

	for (j=1; j<Ny; j++) {

		/***** solve Mij *****/

		// Sij and |Sij|
		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			sr = U.strainrate(i,j,k);
			sra = sqrt( 2.0 * (
					sr[0]*sr[0] + sr[1]*sr[1] + sr[2]*sr[2]
				+ 2.0 * ( sr[3]*sr[3] + sr[4]*sr[4] + sr[5]*sr[5] )
			) );

			S11.id(i,0,k) = sr[0];
			S22.id(i,0,k) = sr[1];
			S33.id(i,0,k) = sr[2];
			S12.id(i,0,k) = sr[3];
			S23.id(i,0,k) = sr[4];
			S13.id(i,0,k) = sr[5];

			EV.id(i,j,k) = sra;
		}}
		// F(Sij), stored in Mij
		fsm11.layerTriFlt(0,0); fsm22.layerTriFlt(0,0); fsm33.layerTriFlt(0,0);
		fsm12.layerTriFlt(0,0); fsm23.layerTriFlt(0,0); fsm13.layerTriFlt(0,0);

		// |Sij|*Sij
		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			sra = EV.id(i,j,k);
			S11.id(i,0,k) *= sra;
			S22.id(i,0,k) *= sra;
			S33.id(i,0,k) *= sra;
			S12.id(i,0,k) *= sra;
			S23.id(i,0,k) *= sra;
			S13.id(i,0,k) *= sra;
		}}
		// F(|Sij|*Sij), stored in Lij
		fsl11.layerTriFlt(0,0); fsl22.layerTriFlt(0,0); fsl33.layerTriFlt(0,0);
		fsl12.layerTriFlt(0,0); fsl23.layerTriFlt(0,0); fsl13.layerTriFlt(0,0);

		dlt1 = pow(   dvol[j], 1./3.);
		dlt2 = pow(4.*dvol[j], 1./3.);

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			// |F(Sij)|
			sra = sqrt(
				M11.id(i,0,k) * M11.id(i,0,k) * 2.
			+	M22.id(i,0,k) * M22.id(i,0,k) * 2.
			+	M33.id(i,0,k) * M33.id(i,0,k) * 2.
			+	M12.id(i,0,k) * M12.id(i,0,k) * 4.
			+	M23.id(i,0,k) * M23.id(i,0,k) * 4.
			+	M13.id(i,0,k) * M13.id(i,0,k) * 4. );
			// Mij
			M11.id(i,0,k) = dlt2*dlt2 * sra * M11.id(i,0,k) - dlt1*dlt1 * L11.id(i,0,k);
			M22.id(i,0,k) = dlt2*dlt2 * sra * M22.id(i,0,k) - dlt1*dlt1 * L22.id(i,0,k);
			M33.id(i,0,k) = dlt2*dlt2 * sra * M33.id(i,0,k) - dlt1*dlt1 * L33.id(i,0,k);
			M12.id(i,0,k) = dlt2*dlt2 * sra * M12.id(i,0,k) - dlt1*dlt1 * L12.id(i,0,k);
			M23.id(i,0,k) = dlt2*dlt2 * sra * M23.id(i,0,k) - dlt1*dlt1 * L23.id(i,0,k);
			M13.id(i,0,k) = dlt2*dlt2 * sra * M13.id(i,0,k) - dlt1*dlt1 * L13.id(i,0,k);
		}}

		/***** solve Lij *****/

		// U in cell centers, temporarily stored in Lii
		U.com1.layerUG2CC(UC1,0,j); // L11
		U.com2.layerVG2CC(UC2,0,j); // L22
		U.com3.layerWG2CC(UC3,0,j); // L33

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			u = UC1.id(i,0,k);
			v = UC2.id(i,0,k);
			w = UC3.id(i,0,k);
			// Ui*Uj
			S11.id(i,0,k) = u * u;
			S22.id(i,0,k) = v * v;
			S33.id(i,0,k) = w * w;
			S12.id(i,0,k) = u * v;
			S23.id(i,0,k) = v * w;
			S13.id(i,0,k) = u * w;
		}}
		// F(Ui*Uj), assigned to Lij, meanwhile filter U and store F(U) in Sij
		fsl12.layerTriFlt(0,0); fsl23.layerTriFlt(0,0); fsl13.layerTriFlt(0,0); // S12->L12, S23->L23, S13->L13
		 flt1.layerTriFlt(0,0);  flt2.layerTriFlt(0,0);  flt3.layerTriFlt(0,0); // L11->S12, L22->S23, L33->S13
		fsl11.layerTriFlt(0,0); fsl22.layerTriFlt(0,0); fsl33.layerTriFlt(0,0); // S11->L11, S22->L22, S33->L33

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			u = UF1.id(i,0,k);
			v = UF2.id(i,0,k);
			w = UF3.id(i,0,k);
			// Lij
			L11.id(i,0,k) = u * u - L11.id(i,0,k);
			L22.id(i,0,k) = v * v - L22.id(i,0,k);
			L33.id(i,0,k) = w * w - L33.id(i,0,k);
			L12.id(i,0,k) = u * v - L12.id(i,0,k);
			L23.id(i,0,k) = v * w - L23.id(i,0,k);
			L13.id(i,0,k) = u * w - L13.id(i,0,k);
			// Lij^d
			sra = 1./3. * (L11.id(i,0,k) + L22.id(i,0,k) + L33.id(i,0,k));
			L11.id(i,0,k) -= sra;
			L22.id(i,0,k) -= sra;
			L33.id(i,0,k) -= sra;
		}}

		/***** calculate eddy viscosity *****/
		dlt1 = dlt2 = 0;

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {

			dlt1 +=	L11.id(i,0,k) * M11.id(i,0,k)
				+	L22.id(i,0,k) * M22.id(i,0,k)
				+	L33.id(i,0,k) * M33.id(i,0,k)
				+	L12.id(i,0,k) * M12.id(i,0,k) * 2.
				+	L23.id(i,0,k) * M23.id(i,0,k) * 2.
				+	L13.id(i,0,k) * M13.id(i,0,k) * 2.;

			dlt2 +=	M11.id(i,0,k) * M11.id(i,0,k)
				+	M22.id(i,0,k) * M22.id(i,0,k)
				+	M33.id(i,0,k) * M33.id(i,0,k)
				+	M12.id(i,0,k) * M12.id(i,0,k) * 2.
				+	M23.id(i,0,k) * M23.id(i,0,k) * 2.
				+	M13.id(i,0,k) * M13.id(i,0,k) * 2.;
		}}

		EV.layerMlt(fmax(0.5 * dlt1 / dlt2, 0) * pow(dvol[j], 2./3.), j);

		// // record for debug
		// FILE *fp = fopen("Cd.dat", "a");
		// char str[64];
		// sprintf(str, "%.18e\n", fmax(0.5 * dlt1 / dlt2, 0));
		// fputs(str, fp);
		// fclose(fp);
	}
	
	EV.layerCpy(EV, 0 , 1);
	EV.layerCpy(EV, Ny, Ny-1);

	mesh0.freeall();
}



# ifdef VREMAN
void SGS::dynamicvreman(Scla &EV, Vctr &U)
{

	int i, j, k;

	double *gr, gr2, beta[6], dy2, Beta;

	for (j=1; j<Ny; j++) {

		dy2 = dy[j] * dy[j];

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {

			gr = U.gradient(i,j,k);
			gr2 =	gr[0]*gr[0] + gr[1]*gr[1] + gr[2]*gr[2]
				+	gr[3]*gr[3] + gr[4]*gr[4] + gr[5]*gr[5]
				+	gr[6]*gr[6] + gr[7]*gr[7] + gr[8]*gr[8];

			beta[0] = dx2 * gr[0]*gr[0] + dy2 * gr[3]*gr[3] + dz2 * gr[6]*gr[6];
			beta[1] = dx2 * gr[1]*gr[1] + dy2 * gr[4]*gr[4] + dz2 * gr[7]*gr[7];
			beta[2] = dx2 * gr[2]*gr[2] + dy2 * gr[5]*gr[5] + dz2 * gr[8]*gr[8];
			beta[3] = dx2 * gr[0]*gr[1] + dy2 * gr[3]*gr[4] + dz2 * gr[6]*gr[7];
			beta[4] = dx2 * gr[1]*gr[2] + dy2 * gr[4]*gr[5] + dz2 * gr[7]*gr[8];
			beta[5] = dx2 * gr[0]*gr[2] + dy2 * gr[3]*gr[5] + dz2 * gr[6]*gr[8];

			Beta =	beta[0]*beta[1] + beta[1]*beta[2] + beta[0]*beta[2]
				-	beta[3]*beta[3] - beta[4]*beta[4] - beta[5]*beta[5];

			EV.id(i,j,k) = Cs * sqrt(Beta / gr2);

		}}
	}
}
# endif




		// // F(Sij), temporarily stored in Lij
		// fsl11.layerTriFlt(0,0); fsl22.layerTriFlt(0,0); fsl33.layerTriFlt(0,0);
		// fsl12.layerTriFlt(0,0); fsl23.layerTriFlt(0,0); fsl13.layerTriFlt(0,0);

		// // |Sij|*Sij
		// S11.layersMlt(q); S22.layersMlt(q); S33.layersMlt(q);
		// S12.layersMlt(q); S23.layersMlt(q); S13.layersMlt(q);

		// // F(|Sij|*Sij), assigned to Mij
		// fsm11.layerTriFlt(0,0); fsm22.layerTriFlt(0,0); fsm33.layerTriFlt(0,0);
		// fsm12.layerTriFlt(0,0); fsm23.layerTriFlt(0,0); fsm13.layerTriFlt(0,0);

		// // - delta_1^2 * F(|Sij|*Sij)
		// dlt = - pow(dvol[j], 2./3.);
		// M11.layerMlt(dlt); M22.layerMlt(dlt); M33.layerMlt(dlt);
		// M12.layerMlt(dlt); M23.layerMlt(dlt); M13.layerMlt(dlt);

		// // |F(Sij)|
		// for (k=0; k<Nz; k++) {
		// for (i=0; i<Nx; i++) {
		// 	q.id(i,0,k) = sqrt(
		// 		L11.id(i,0,k) * L11.id(i,0,k) * 2.
		// 	+	L22.id(i,0,k) * L22.id(i,0,k) * 2.
		// 	+	L33.id(i,0,k) * L33.id(i,0,k) * 2.
		// 	+	L12.id(i,0,k) * L12.id(i,0,k) * 4.
		// 	+	L23.id(i,0,k) * L23.id(i,0,k) * 4.
		// 	+	L13.id(i,0,k) * L13.id(i,0,k) * 4. );
		// }}

		// // delta_2^2 * |F(Sij)|*F(Sij)
		// dlt = pow(4.*dvol[j], 2./3.);
		// L11.layersMlt(q).layerMlt(dlt);
		// L22.layersMlt(q).layerMlt(dlt);
		// L33.layersMlt(q).layerMlt(dlt);
		// L12.layersMlt(q).layerMlt(dlt);
		// L23.layersMlt(q).layerMlt(dlt);
		// L13.layersMlt(q).layerMlt(dlt);

		// // Mij = delta_2^2 * |F(Sij)|*F(Sij) - delta_1^2 * F(|Sij|*Sij)
		// M11.layersAdd(L11); M22.layersAdd(L22); M33.layersAdd(L33);
		// M12.layersAdd(L12); M23.layersAdd(L23); M13.layersAdd(L13);

		
		// // Ui*Uj, temporarily stored in Sij
		// S11.layerCpy(UC1).layersMlt(UC1);
		// S22.layerCpy(UC2).layersMlt(UC2);
		// S33.layerCpy(UC3).layersMlt(UC3);
		// S12.layerCpy(UC1).layersMlt(UC2);
		// S23.layerCpy(UC2).layersMlt(UC3);
		// S13.layerCpy(UC1).layersMlt(UC3);

		// // Lij = Ui*Uj - F(Ui*Uj)
		// L11.layerMlt(-1); L22.layerMlt(-1); L33.layerMlt(-1);
		// L12.layerMlt(-1); L23.layerMlt(-1); L13.layerMlt(-1);

		// L11.layersAdd(q.layerCpy(UF1).layersMlt(UF1));
		// L22.layersAdd(q.layerCpy(UF2).layersMlt(UF2));
		// L33.layersAdd(q.layerCpy(UF3).layersMlt(UF3));
		// L12.layersAdd(q.layerCpy(UF1).layersMlt(UF2));
		// L23.layersAdd(q.layerCpy(UF2).layersMlt(UF3));
		// L13.layersAdd(q.layerCpy(UF1).layersMlt(UF3));

		// // Lij^d = Lij - 1/3 * delta_ij * Lkk
		// q.layerCpy(L11).layersAdd(L22).layersAdd(L33).layerMlt(-1./3.);
		// L11.layersAdd(q);
		// L22.layersAdd(q);
		// L33.layersAdd(q);



# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>

# include "Interp.h"
# include "SGS.h"


// # define SGSDEBUG



void SGS::smargorinsky(double Cs, double Re)
{
	int i, j, k;
	double *sr, sra, dlt, dmp;

	// calculate delta_nu. NOTE: this only applys when boundary is on the wall
	double tauw = (0.5/Re) * (
		( U[1].av(1) - U[1].av(0)    ) / h[1]
	  - ( U[1].av(Ny)- U[1].av(Ny-1) ) / h[Ny] );
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
	EV.lyrSet(EV[1], 0).lyrSet(EV[Ny-1], Ny);
}


void SGS::dynamicsmarg()
{
	int i, j, k;
	double *sr, sra, dlt1, dlt2, u, v, w;

# ifdef SGSDEBUG
	FILE *fp = fopen("Cd.dat", "w");
	char str[64];
# endif

	Mesh ms(Nx,0,Nz,Lx,0,Lz);
	Scla S11(ms), S22(ms), S33(ms), S12(ms), S23(ms), S13(ms);
	Scla M11(ms), M22(ms), M33(ms), M12(ms), M23(ms), M13(ms);
	Scla L11(ms), L22(ms), L33(ms), L12(ms), L23(ms), L13(ms);
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

			S11.id(i,0,k) = sr[0]; S22.id(i,0,k) = sr[1]; S33.id(i,0,k) = sr[2];
			S12.id(i,0,k) = sr[3]; S23.id(i,0,k) = sr[4]; S13.id(i,0,k) = sr[5];

			EV.id(i,j,k) = sra;
		}}
		// F(Sij), stored in Mij
		fsm11.layerTriFlt(0,0); fsm22.layerTriFlt(0,0); fsm33.layerTriFlt(0,0);
		fsm12.layerTriFlt(0,0); fsm23.layerTriFlt(0,0); fsm13.layerTriFlt(0,0);

		// |Sij|*Sij
		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			sra = EV.id(i,j,k);
			S11.id(i,0,k) *= sra; S22.id(i,0,k) *= sra; S33.id(i,0,k) *= sra;
			S12.id(i,0,k) *= sra; S23.id(i,0,k) *= sra; S13.id(i,0,k) *= sra;
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
		U.layer2CC(UC1[0], UC2[0], UC3[0], j); // L11, L22, L33
		// U.com1.layerUG2CC(UC1,0,j); // L11
		// U.com2.layerVG2CC(UC2,0,j); // L22
		// U.com3.layerWG2CC(UC3,0,j); // L33

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			u = UC1.id(i,0,k);
			v = UC2.id(i,0,k);
			w = UC3.id(i,0,k);
			// Ui*Uj
			S11.id(i,0,k) = u * u; S22.id(i,0,k) = v * v; S33.id(i,0,k) = w * w;
			S12.id(i,0,k) = u * v; S23.id(i,0,k) = v * w; S13.id(i,0,k) = u * w;
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
		for (i=0; i<Nx; i++) {
			dlt1 = dlt2 = 0;

			for (k=0; k<Nz; k++) {

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
			}

			for (k=0; k<Nz; k++)
				EV.id(i,j,k) *= fmin(fmax(.5*dlt1/dlt2, 0), .5) * pow(dvol[j], 2./3.); // |Sij| has been assigned to EV

# ifdef SGSDEBUG
			sprintf(str, "%.18e\t", fmin(fmax(.5*dlt1/dlt2, 0), .5));
			fputs(str, fp);
# endif
			
		}

# ifdef SGSDEBUG
		fputs("\n", fp);
		if (j==Ny-1) fclose(fp);
# endif

	}
	
	EV.lyrSet(EV[1], 0).lyrSet(EV[Ny-1], Ny);
	ms.freeall();
}



void SGS::dynamicvreman(double Re)
{
	int i, j, k;
	double *sr, fsr[6], sr2;
	double *gr, fgr[9], gr2;
	double dy2, beta[6], Beta;
	double sum1 = 0, sum2 = 0;

# ifdef SGSDEBUG
	FILE *fp = fopen("Cv.dat", "a");
	char str[64];
# endif

	Mesh ms(Nx,0,Nz,Lx,0,Lz);

	Scla	FG11(ms), FG12(ms), FG13(ms), \
			FG21(ms), FG22(ms), FG23(ms), \
			FG31(ms), FG32(ms), FG33(ms), \
			FGG (ms), FPSS(ms), FSFS(ms);

	Scla	G11(ms), G12(ms), G13(ms), \
			G21(ms), G22(ms), G23(ms), \
			G31(ms), G32(ms), G33(ms), \
			GG (ms), PSS(ms), &SS = FSFS, \
			&S11 = G11, &S22 = G12, &S33 = G13, \
			&S12 = G21, &S23 = G22, &S13 = G23, \
			&FS11 = FG11, &FS22 = FG12, &FS33 = FG13, \
			&FS12 = FG21, &FS23 = FG22, &FS13 = FG23;

	Interp	fg11(G11,FG11), fg12(G12,FG12), fg13(G13,FG13), \
			fg21(G21,FG21), fg22(G22,FG22), fg23(G23,FG23), \
			fg31(G31,FG31), fg32(G32,FG32), fg33(G33,FG33), \
			fpss(PSS,FPSS), fgg(GG,FGG), \
			&fs11 = fg11, &fs22 = fg12, &fs33 = fg13, \
			&fs12 = fg21, &fs23 = fg22, &fs13 = fg23;


	for (j=1; j<Ny; j++) {
		dy2 = dy[j] * dy[j];
		// Sij, SS
		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			sr = U.strainrate(i,j,k);
			sr2 =	sr[0]*sr[0] + sr[1]*sr[1] + sr[2]*sr[2]
				+ (	sr[3]*sr[3] + sr[4]*sr[4] + sr[5]*sr[5] ) * 2.;
			S11.id(i,0,k) = sr[0]; S22.id(i,0,k) = sr[1]; S33.id(i,0,k) = sr[2]; // stroed in FG1i
			S12.id(i,0,k) = sr[3]; S23.id(i,0,k) = sr[4]; S13.id(i,0,k) = sr[5]; // stored in FG2i
			SS.id(i,0,k) = sr2;
		}}
		// F(Sij)
		fs11.layerTriFlt(0,0); fs22.layerTriFlt(0,0); fs33.layerTriFlt(0,0); // FG1i -> FSii
		fs12.layerTriFlt(0,0); fs23.layerTriFlt(0,0); fs13.layerTriFlt(0,0); // FG2i -> FSij

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			fsr[0] = FS11.id(i,0,k); fsr[1] = FS22.id(i,0,k); fsr[2] = FS33.id(i,0,k);
			fsr[3] = FS12.id(i,0,k); fsr[4] = FS23.id(i,0,k); fsr[5] = FS13.id(i,0,k);
			// dUj/dXi
			gr = U.gradient(i,j,k);
			G11.id(i,0,k) = gr[0]; G12.id(i,0,k) = gr[1]; G13.id(i,0,k) = gr[2];
			G21.id(i,0,k) = gr[3]; G22.id(i,0,k) = gr[4]; G23.id(i,0,k) = gr[5];
			G31.id(i,0,k) = gr[6]; G32.id(i,0,k) = gr[7]; G33.id(i,0,k) = gr[8];

			sr2 =	fsr[0]*fsr[0] + fsr[1]*fsr[1] + fsr[2]*fsr[2]
				+ (	fsr[3]*fsr[3] + fsr[4]*fsr[4] + fsr[5]*fsr[5] ) * 2.;
			
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
				- (	beta[3]*beta[3] + beta[4]*beta[4] + beta[5]*beta[5]	);
	
			// PI^g, PI^g * SS, (dUj/dXi)(dUj/dXi), F(Sij)F(Sij)
			PSS.id(i,0,k) = SS.id(i,0,k) * ( EV.id(i,j,k) = sqrt(Beta / gr2) );
			GG.id(i,0,k) = gr2;
			FSFS.id(i,0,k) = sr2;
		}}
		// F(dUj/dXi)
		fg11.layerTriFlt(0,0); fg12.layerTriFlt(0,0); fg13.layerTriFlt(0,0);
		fg21.layerTriFlt(0,0); fg22.layerTriFlt(0,0); fg23.layerTriFlt(0,0);
		fg31.layerTriFlt(0,0); fg32.layerTriFlt(0,0); fg33.layerTriFlt(0,0);
		// F( (dUj/dXi)(dUj/dXi) )
		fgg.layerTriFlt(0,0);
		// F( PI^g * SijSij )
		fpss.layerTriFlt(0,0);
	
		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			// F(dUj/dXi)F(dUj/dXi)
			fgr[0] = FG11.id(i,0,k); fgr[1] = FG12.id(i,0,k); fgr[2] = FG13.id(i,0,k);
			fgr[3] = FG21.id(i,0,k); fgr[4] = FG22.id(i,0,k); fgr[5] = FG23.id(i,0,k);
			fgr[6] = FG31.id(i,0,k); fgr[7] = FG32.id(i,0,k); fgr[8] = FG33.id(i,0,k);
			gr2 =	fgr[0]*fgr[0] + fgr[1]*fgr[1] + fgr[2]*fgr[2]
				+	fgr[3]*fgr[3] + fgr[4]*fgr[4] + fgr[5]*fgr[5]
				+	fgr[6]*fgr[6] + fgr[7]*fgr[7] + fgr[8]*fgr[8];
			// PI^t
			beta[0] = 4.*dx2 * fgr[0]*fgr[0] + dy2 * fgr[3]*fgr[3] + 4.*dz2 * fgr[6]*fgr[6];
			beta[1] = 4.*dx2 * fgr[1]*fgr[1] + dy2 * fgr[4]*fgr[4] + 4.*dz2 * fgr[7]*fgr[7];
			beta[2] = 4.*dx2 * fgr[2]*fgr[2] + dy2 * fgr[5]*fgr[5] + 4.*dz2 * fgr[8]*fgr[8];
			beta[3] = 4.*dx2 * fgr[0]*fgr[1] + dy2 * fgr[3]*fgr[4] + 4.*dz2 * fgr[6]*fgr[7];
			beta[4] = 4.*dx2 * fgr[1]*fgr[2] + dy2 * fgr[4]*fgr[5] + 4.*dz2 * fgr[7]*fgr[8];
			beta[5] = 4.*dx2 * fgr[0]*fgr[2] + dy2 * fgr[3]*fgr[5] + 4.*dz2 * fgr[6]*fgr[8];

			Beta =	beta[0]*beta[1] + beta[1]*beta[2] + beta[0]*beta[2]
				- (	beta[3]*beta[3] + beta[4]*beta[4] + beta[5]*beta[5]	);

			sum1 += dvol[j] * ( FGG.id(i,0,k) - gr2 );
			sum2 += dvol[j] * ( FPSS.id(i,0,k) - sqrt(Beta / gr2) * FSFS.id(i,0,k) );
		}}
	}

	EV *= fmin(fmax(-.5/Re*sum1/sum2, 0), .5); // PI^g has been assigned to EV
	// EV.bulkMlt(fmin(fmax(-.5/Re*sum1/sum2, 0), .5));
	// EV.bulkMlt(0.05);

	EV.lyrSet(EV[1], 0).lyrSet(EV[Ny-1], Ny);
	ms.freeall();

# ifdef SGSDEBUG
	sprintf(str, "%.18e\n", fmin(fmax(-.5/Re*sum1/sum2, 0), .5));
	fputs(str, fp);
	fclose(fp);
# endif

}



// double SGS::vreman(Vctr &U, int i, int j, int k)
// {
// 	double *gr, gr2, beta[6], Beta, dy2 = dy[j]*dy[j];

// 	gr = U.gradient(i,j,k);
// 	gr2 =	gr[0]*gr[0] + gr[1]*gr[1] + gr[2]*gr[2]
// 		+	gr[3]*gr[3] + gr[4]*gr[4] + gr[5]*gr[5]
// 		+	gr[6]*gr[6] + gr[7]*gr[7] + gr[8]*gr[8];

// 	beta[0] = dx2 * gr[0]*gr[0] + dy2 * gr[3]*gr[3] + dz2 * gr[6]*gr[6];
// 	beta[1] = dx2 * gr[1]*gr[1] + dy2 * gr[4]*gr[4] + dz2 * gr[7]*gr[7];
// 	beta[2] = dx2 * gr[2]*gr[2] + dy2 * gr[5]*gr[5] + dz2 * gr[8]*gr[8];
// 	beta[3] = dx2 * gr[0]*gr[1] + dy2 * gr[3]*gr[4] + dz2 * gr[6]*gr[7];
// 	beta[4] = dx2 * gr[1]*gr[2] + dy2 * gr[4]*gr[5] + dz2 * gr[7]*gr[8];
// 	beta[5] = dx2 * gr[0]*gr[2] + dy2 * gr[3]*gr[5] + dz2 * gr[6]*gr[8];

// 	Beta =	beta[0]*beta[1] + beta[1]*beta[2] + beta[0]*beta[2]
// 		- (	beta[3]*beta[3] + beta[4]*beta[4] + beta[5]*beta[5]	);

// 	return sqrt(Beta / gr2);
// }









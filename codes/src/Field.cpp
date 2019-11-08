# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <cmath>

# include "Field.h"
# include "Interp.h"

using namespace std;


/***** fields initiation *****/

void Field::initField(double energy)
/* initiate flow field from laminar with random fluctions */
{
	int idx, i, j, k;
	Scla &U1 = U.com1, &U2 = U.com2, &U3 = U.com3;

	U1.bulkSet(0);
	U2.bulkSet(0);
	U3.bulkSet(0);
	 P.bulkSet(0);

	// initiate velocities with random fluctuations
	srand(time(0));
	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		U1.id(i,j,k) = energy * (rand()-rand()) / (double)(RAND_MAX); // probability distribution: p(x) = ( x<0 ? x+1 : 1-x )
		if (j>1)
		U2.id(i,j,k) = energy * (rand()-rand()) / (double)(RAND_MAX); // 2nd order moment is 2/3, <U^2+V^2+W^2>/2 = energy
		U3.id(i,j,k) = energy * (rand()-rand()) / (double)(RAND_MAX); // sample space expands on the whole physical domain
	}}}

	// remove mean velocities
	for (j=2; j<Ny; j++) {	// remove mean V at every XZ plane
		U2.layerAdd(- U2.layerMean(j), j);
	}
	for (i=0; i<Nx; i++) {	// remove mean U at every YZ plane
		double mean = 0.0;
		for (k=0; k<Nz; k++) mean += U1.yMeanU(i, k) / Nz;
		for (j=1; j<Ny; j++) {
		for (k=0; k<Nz; k++) {
			U1.id(i,j,k) -= mean;
		}}
	}
	for (k=0; k<Nz; k++) {	// remove mean W at every XY plane
		double mean = 0.0;
		for (i=0; i<Nx; i++) mean += U3.yMeanU(i, k) / Nx;
		for (j=1; j<Ny; j++) {
		for (i=0; i<Nx; i++) {
			U3.id(i,j,k) -= mean;
		}}
	}

	// impose parabolic profile from laminar flow (Ly is taken for 2.0)
	for (j=1; j<Ny; j++) U1.layerAdd(yc[j] * (2.-yc[j]), j);

	// modify flow rate // note: boundary will be modified through using bulk functions
	U1.bulkMlt(1.0 / U1.bulkMeanU()); // rescale the mass flow rate to be equal to 2.0 (bulk mean U = 1.0 because of non-dimensionalization)
	U3.bulkAdd(- U3.bulkMeanU());

	// implement BC on boundaries
	this->bcond(0);
	this->applyBC();
}

void Field::initField(Field &field)
{
	Scla *src[4] = {&field.U.com1, &field.U.com2, &field.U.com3, &field.P};
	Scla *dst[4] = {&U.com1, &U.com2, &U.com3, &P};

	for (int n=0; n<4; n++) Interp(*src[n], *dst[n]).interpolate(n==1?'V':'U');
	
	// for (int n=0; n<4; n++) {
	// 	field.DP.bulkCpy(*src[n]);
	// 	dst[n]->bulkCpy(DP.interpolate(field.DP));
	// }
}


/***** boundary (UBC) process *****/

void Field::bcond(int tstep)
/* set boundary conditions (not yet applied to velocity field) */
{
	UBC.com1.bulkSet(0);
	UBC.com2.bulkSet(0);
	UBC.com3.bulkSet(0);
}

void Field::bcond(Vctr &UBC0)
{
	UBC.com1.bulkCpy(UBC0.com1);
	UBC.com2.bulkCpy(UBC0.com2);
	UBC.com3.bulkCpy(UBC0.com3);
}


/***** computation related functions *****/

void Field::getup(double dt, int nthrds)
{
	this->idm.ompset(nthrds);
	this->idm.dtset(dt);
	
	this->idm.uhcalc(UH, U , P, UBC, NU, mpg);
	this->idm.dpcalc(DP, UH, P, UBC         );
	this->idm.upcalc(U , P , UH, DP, mpg    );
	this->applyBC(dt);
}

void Field::getnu(double Re, double Cs)
{
	if (Cs > 0) {
		this->sgs.evcalc(NU, U, Re, Cs);
		NU.bulkAdd(1.0 / Re);
	}
	else  NU.bulkSet(1.0 / Re);
}

void Field::applyBC()
/* apply Dirichlet BC on velocities */
{
	U.com1.layerCpy(UBC.com1, 0, 0);
	U.com2.layerCpy(UBC.com2, 1, 0);
	U.com3.layerCpy(UBC.com3, 0, 0);
	U.com1.layerCpy(UBC.com1, Ny, 1);
	U.com2.layerCpy(UBC.com2, Ny, 1);
	U.com3.layerCpy(UBC.com3, Ny, 1);
}

void Field::applyBC(double dt)
/* apply Dirichlet BC on velocities, with boundary time derivative considered */
{
	Scla &U1 = U.com1,   &U2 = U.com2,   &U3 = U.com3;
	Scla &H1 = UH.com1,  &H2 = UH.com2,  &H3 = UH.com3;
	Scla &B1 = UBC.com1, &B2 = UBC.com2, &B3 = UBC.com3;
	Scla ul(Mesh(Nx,0,Nz,Lx,0,Lz)), vl(ul.meshGet()), wl(ul.meshGet());
	
	ul.layerCpy(U1,0,0).layerMlt(-1); ul.layersAdd(B1,0,0).layerMlt(1.0/dt);
	vl.layerCpy(U2,0,1).layerMlt(-1); vl.layersAdd(B2,0,0).layerMlt(1.0/dt);
	wl.layerCpy(U3,0,0).layerMlt(-1); wl.layersAdd(B3,0,0).layerMlt(1.0/dt);
	H1.layerCpy(ul, 0);
	H2.layerCpy(vl, 1);
	H3.layerCpy(wl, 0);

	ul.layerCpy(U1,0,Ny).layerMlt(-1); ul.layersAdd(B1,0,1).layerMlt(1.0/dt);
	vl.layerCpy(U2,0,Ny).layerMlt(-1); vl.layersAdd(B2,0,1).layerMlt(1.0/dt);
	wl.layerCpy(U3,0,Ny).layerMlt(-1); wl.layersAdd(B3,0,1).layerMlt(1.0/dt);
	H1.layerCpy(ul, Ny);
	H2.layerCpy(vl, Ny);
	H3.layerCpy(wl, Ny);

	ul.meshGet().freeall();
	this->applyBC();
}

void Field::removeSpanMean()
{
	int i, j, k, n;
	double qm, *qsm = new double [Nx];
	Bulk *bks[4] = {&U.com1, &U.com2, &UH.com1, &UH.com2};

	for (n=0; n<4; n++) {
	for (j=1; j<Ny; j++) {
		if (j==1 && (n==1 || n==3)) continue; // j==1 layer for V & VH is boundary, skip

		for (i=0; i<Nx; i++) { qsm[i] = 0.0;
		for (k=0; k<Nz; k++) { qsm[i] += bks[n]->id(i,j,k) / Nz; }}

		qm = 0.0;
		for (i=0; i<Nx; i++)   qm += qsm[i] / Nx;

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) { bks[n]->id(i,j,k) -= qsm[i] - qm; }}
	}}

	delete [] qsm;
}


/***** file IO operations *****/

void Field::writeField(char *path, int tstep)
{
	Bulk *bks[8] = {
		&U.com1,  &U.com2,  &U.com3,  &P,
		&UH.com1, &UH.com2, &UH.com3, &DP };
	char names[8][32] = {"U", "V", "W", "P", "UT", "VT", "WT", "PT"};
	for (int n=0; n<8; n++) {
		sprintf(names[n], "%s%08i", names[n], tstep);
		bks[n]->fileIO(path, names[n], 'w');
	}
}

void Field::readField(char *path, int tstep)
{
	Bulk *bks[4] = {&U.com1, &U.com2, &U.com3, &P};
	char names[4][32] = {"U", "V", "W", "P"};
	for (int n=0; n<4; n++) {
		sprintf(names[n], "%s%08i", names[n], tstep);
		bks[n]->fileIO(path, names[n], 'r');
	}
}

void Field::writeTecplot(char *path, int tstep, double time)
/* write velocity and pressure fields to ascii files readable by tecplot */
{
	FILE *fp;
	char str[1024], fn[1024];
	int i, j, k, idx, ip, jp, kp;
	double v1, v2, *u = U.bulkGet(1), *v = U.bulkGet(2), *w = U.bulkGet(3), *p = P.bulkGet();

	sprintf(fn, "%sFIELD%08i.dat", path?path:"", tstep);

	fp = fopen(fn, "w");

		fputs("Title = \"3D instantaneous field\"\n", fp);
		fputs("variables = \"x\", \"z\", \"y\", \"u\", \"v\", \"w\", \"p\"\n", fp);
		sprintf(str, "zone t = \"%f\", i = %i, j = %i, k = %i\n", time, Ny+1, Nz, Nx);
		fputs(str, fp);

		for (i=0; i<Nx; i++) {
		for (k=0; k<Nz; k++) {
		for (j=0; j<=Ny;j++) {
			idx = IDX(i,j,k);
			ip = IDX(ipa[i],j,k);
			jp = IDX(i,j+1,k);
			kp = IDX(i,j,kpa[k]);

			v1 = j==0 ? v[jp] : v[idx];
			v2 = j==Ny ? v[idx] : v[jp];
			// all interpolate to cell centers. tool functions not used here because of the order of i,j,k
			sprintf(str, "%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
				dx * (i + 0.5),
				dz * (k + 0.5),
				yc[j],
				0.5 * (u[idx] + u[ip]),
				0.5 * (v1 + v2),
				0.5 * (w[idx] + w[kp]),
				p[idx]	);
			fputs(str, fp);
		}}}

	fclose(fp);

	sprintf(str, "preplot %s", fn);	system(str);
	sprintf(str, "rm %s", fn);		system(str);
}

void Field::debug_Output(int tstep)
{
	char path[1024] = "debug/";
	{
		char names[4][32] = {"U", "V", "W", "P"};
		U.com1.debug_AsciiOutput(path, names[0], 0, Ny+1);
		U.com2.debug_AsciiOutput(path, names[1], 1, Ny+1);
		U.com3.debug_AsciiOutput(path, names[2], 0, Ny+1);
		     P.debug_AsciiOutput(path, names[3], 1, Ny);
	}{
		char names[4][32] = {"UH", "VH", "WH", "DP"};
		UH.com1.debug_AsciiOutput(path, names[0], 0, Ny+1);
		UH.com2.debug_AsciiOutput(path, names[1], 1, Ny+1);
		UH.com3.debug_AsciiOutput(path, names[2], 0, Ny+1);
		     DP.debug_AsciiOutput(path, names[3], 1, Ny);
	}{
		char names[4][32] = {"UBC", "VBC", "WBC", "NU"};
		UBC.com1.debug_AsciiOutput(path, names[0], 0, 2);
		UBC.com2.debug_AsciiOutput(path, names[1], 1, 2);
		UBC.com3.debug_AsciiOutput(path, names[2], 0, 2);
		      NU.debug_AsciiOutput(path, names[3], 1, Ny);
	}
}







// # define DEBUG	// g++ -lfftw3 -lm -I include src/Field.cpp src/Mesh.cpp
# ifdef DEBUG

int main()
{
	int dim[3] = {8, 2, 12}, Ny = dim[1], Nxz = dim[0] * dim[2];
	class Field *pfield = new class Field(dim);

	double *mat = & (pfield->DP[Nxz]);
	mat[0] = 1.0;
	for (int k=0; k<dim[2]; k++) {
	for (int i=0; i<dim[0]; i++) {
		mat[dim[0]*k+i] = rand() / (double)(RAND_MAX);
	}}


	pfield->fft();
	pfield->debug_AsciiOutputF(0, 1, Ny-1);

	pfield->ifft();
	pfield->debug_AsciiOutput(0, pfield->DP, 1, Ny-1, "DP");
}

# endif








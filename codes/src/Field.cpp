# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <cmath>

# include "Field.h"
# include "Interp.h"

using namespace std;



Field::Field(const Mesh &mesh):
Mesh(mesh),
U(mesh), UH(mesh),   FB(mesh),
P(mesh), DP(mesh,1), NU(mesh),
UBC(Mesh(Nx,2-1,Nz,Lx,Ly,Lz)),
PBC(Mesh(Nx,2-1,Nz,Lx,Ly,Lz))
{
	UBC.meshGet().y[0] = 0;
	UBC.meshGet().y[1] = 0;
	UBC.meshGet().yc[0] = y[1];
	UBC.meshGet().yc[1] = y[Ny];
}


/***** fields initiation *****/

void Field::initField(double energy)
/* initiate flow field (U, V, W, P, including boundaries) from laminar with random fluctions */
{
	int i, j, k;

	U[1] = 0.; U[2] = 0.; U[3] = 0.; P = 0.;

	// initiate velocities with random fluctuations
	srand(time(0));
	for (j=1; j<Ny; j++) { // probability distribution: p(x) = ( x<0 ? x+1 : 1-x )
	for (k=0; k<Nz; k++) { // 2nd order moment is 2/3, <U^2+V^2+W^2>/2 = energy
	for (i=0; i<Nx; i++) { // sample space expands on the whole physical domain
		U[1].id(i,j,k) = energy * (rand()-rand()) / (double)(RAND_MAX); if (j>1) // act on nex tline
		U[2].id(i,j,k) = energy * (rand()-rand()) / (double)(RAND_MAX);
		U[3].id(i,j,k) = energy * (rand()-rand()) / (double)(RAND_MAX);
	}}}

	// remove mean velocities
	for (j=2; j<Ny; j++) {	// remove mean V of every XZ plane
		U[2].lyrAdd(- U[2].av(j), j);
	}
	for (i=0; i<Nx; i++) {	// remove mean U of every YZ plane
		double mean = 0.0;
		for (k=0; k<Nz; k++) mean += U[1].yMeanU(i, k) / Nz;
		for (j=1; j<Ny; j++) {
		for (k=0; k<Nz; k++) {
			U[1].id(i,j,k) -= mean;
		}}
	}
	for (k=0; k<Nz; k++) {	// remove mean W of every XY plane
		double mean = 0.0;
		for (i=0; i<Nx; i++) mean += U[3].yMeanU(i, k) / Nx;
		for (j=1; j<Ny; j++) {
		for (i=0; i<Nx; i++) {
			U[3].id(i,j,k) -= mean;
		}}
	}

	// impose parabolic profile from laminar flow (Ly is taken for 2.0)
	for (j=1; j<Ny; j++) U[1].lyrAdd(yc[j] * (2.-yc[j]), j);

	// modify flow rate
	U[1] *= 1. / U[1].bulkMeanU(); // bulk mean U = 1.0 due to non-dimensionalization
	U[3] += - U[3].bulkMeanU(); // note: boundaries are modified through using bulk functions

	// initiate boundaries with homogeneous BC
	U[1].lyrSet(U[1][1], 0).lyrMlt(- dy[0]/dy[1], 0).lyrSet(U[1][Ny-1], Ny).lyrMlt(- dy[Ny]/dy[Ny-1], Ny);
	U[3].lyrSet(U[3][1], 0).lyrMlt(- dy[0]/dy[1], 0).lyrSet(U[3][Ny-1], Ny).lyrMlt(- dy[Ny]/dy[Ny-1], Ny);
	U[2].lyrSet(0., 1).lyrSet(0., Ny);

	// set all other variables to 0
	this->reset();
}

void Field::initField(Field &field)
/* initiate flow field (U, V, W, P, including boundaries) from existing fields */
{
	Scla *src[4] = {&field.U[1], &field.U[2], &field.U[3], &field.P};
	Scla *dst[4] = {&      U[1], &      U[2], &      U[3], &      P};
	for (int n=0; n<4; n++) Interp(*src[n], *dst[n]).bulkInterp(n==1?'V':'U');
	this->reset();
}

void Field::reset()
{
	U[2].lyrSet(0., 0); // j=0 layer of V does not participate in computation
	DP = 0.; NU = 0.; PBC = 0.;
	UH [1] = 0.; UH [2] = 0.; UH [3] = 0.;
	FB [1] = 0.; FB [2] = 0.; FB [3] = 0.;
	UBC[1] = 0.; UBC[2] = 0.; UBC[3] = 0.;
	mpg[0] = 0.; mpg[1] = 0.; mpg[2] = 0.;
}


/***** file IO operations *****/

void Field::writeField(char *path, int tstep)
{
	Bulk *bks[8] = {
		&U[1],  &U[2],  &U[3],  &P,
		&UH[1], &UH[2], &UH[3], &DP };
	char names[8][32] = {"U", "V", "W", "P", "UT", "VT", "WT", "PT"};
	for (int n=0; n<8; n++) {
		sprintf(names[n], "%s%08i", names[n], tstep);
		bks[n]->fileIO(path, names[n], 'w');
	}
}

Field& Field::readField(char *path, int tstep)
{
	Bulk *bks[4] = {&U[1], &U[2], &U[3], &P};
	char names[4][32] = {"U", "V", "W", "P"};
	for (int n=0; n<4; n++) {
		sprintf(names[n], "%s%08i", names[n], tstep);
		bks[n]->fileIO(path, names[n], 'r');
	}
	return (Field&)(*this);
}

void Field::writeTecplot(char *path, int tstep, double time)
/* write velocity and pressure fields to ascii files readable by tecplot */
{
	FILE *fp;
	char str[1024], fn[1024];
	int i, j, k;
	// int idx, ip, jp, kp;
	// double v1, v2, *u = U[1].blkGet(), *v = U[2].blkGet(), *w = U[3].blkGet(), *p = P.blkGet();

	sprintf(fn, "%sFIELD%08i.dat", path?path:"", tstep);

	fp = fopen(fn, "w");


	fputs("Title = \"3D instantaneous field\"\n", fp);
	fputs("variables = \"x\", \"y\", \"z\", \"u\", \"v\", \"w\", \"p\"\n", fp);
	sprintf(str, "zone t = \"%f\", i = %i, j = %i, k = %i\n", time, Nx, Nz, Ny+1);
	fputs(str, fp);

	Mesh ms(Nx,0,Nz,Lx,0,Lz);
	Scla ul(ms), vl(ms), wl(ms), pl(ms);

	for (j=0; j<=Ny; j++) {

		U.layer2CC(ul[0], vl[0], wl[0], j);
		pl = P[j];

		for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			sprintf(str, "%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
				dx * (i + 0.5),
				j==0 ? y[1] : j==Ny ? y[Ny] : yc[j],
				dz * (k + 0.5),
				ul.id(i,0,k),
				vl.id(i,0,k),
				wl.id(i,0,k),
				pl.id(i,0,k)
			);
			fputs(str, fp);
		}}
	}

	ms.freeall();


		// fputs("Title = \"3D instantaneous field\"\n", fp);
		// fputs("variables = \"x\", \"z\", \"y\", \"u\", \"v\", \"w\", \"p\"\n", fp);
		// sprintf(str, "zone t = \"%f\", i = %i, j = %i, k = %i\n", time, Ny+1, Nz, Nx);
		// fputs(str, fp);

		// for (i=0; i<Nx; i++) {
		// for (k=0; k<Nz; k++) {
		// for (j=0; j<=Ny;j++) {
		// 	idx = IDX(i,j,k);
		// 	ip = IDX(ipa[i],j,k);
		// 	jp = IDX(i,j+1,k);
		// 	kp = IDX(i,j,kpa[k]);

		// 	v1 = j==0 ? v[jp] : v[idx];
		// 	v2 = j==Ny ? v[idx] : v[jp];
		// 	// all interpolate to cell centers. tool functions not used here because of the order of i,j,k
		// 	sprintf(str, "%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
		// 		dx * (i + 0.5),
		// 		dz * (k + 0.5),
		// 		yc[j],
		// 		0.5 * (u[idx] + u[ip]),
		// 		0.5 * (v1 + v2),
		// 		0.5 * (w[idx] + w[kp]),
		// 		p[idx]	);
		// 	fputs(str, fp);
		// }}}

	fclose(fp);

	sprintf(str, "preplot %s", fn);	system(str);
	sprintf(str, "rm %s", fn);		system(str);
}

void Field::debug_Output(int tstep)
{
	char path[1024] = "debug/";
	{
		char names[4][32] = {"U", "V", "W", "P"};
		U[1].debug_AsciiOutput(path, names[0], 0, Ny+1);
		U[2].debug_AsciiOutput(path, names[1], 1, Ny+1);
		U[3].debug_AsciiOutput(path, names[2], 0, Ny+1);
		   P.debug_AsciiOutput(path, names[3], 1, Ny);
	}{
		char names[4][32] = {"UH", "VH", "WH", "DP"};
		UH[1].debug_AsciiOutput(path, names[0], 0, Ny+1);
		UH[2].debug_AsciiOutput(path, names[1], 1, Ny+1);
		UH[3].debug_AsciiOutput(path, names[2], 0, Ny+1);
		   DP.debug_AsciiOutput(path, names[3], 1, Ny);
	}{
		char names[4][32] = {"UBC", "VBC", "WBC", "NU"};
		UBC[1].debug_AsciiOutput(path, names[0], 0, 2);
		UBC[2].debug_AsciiOutput(path, names[1], 1, 2);
		UBC[3].debug_AsciiOutput(path, names[2], 0, 2);
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








# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>
# include <time.h>

# include "Basic.h"
# include "Interp.h"

using namespace std;


Feld& Feld::initrand(double energy)
/* initiate flow field (U, V, W, P, including boundaries) from laminar with random fluctions */
{
	int i, j, k;

	this->reset(0);

	// initiate velocities with random fluctuations
	srand(time(0));
	for (j=1; j<Ny; j++) { // probability distribution: p(x) = ( x<0 ? x+1 : 1-x )
	for (k=0; k<Nz; k++) { // 2nd order moment is 2/3, <U^2+V^2+W^2>/2 = energy
	for (i=0; i<Nx; i++) { // sample space expands on the whole physical domain
		V[1].id(i,j,k) = energy * (rand()-rand()) / (double)(RAND_MAX); if (j>1) // act on nex tline
		V[2].id(i,j,k) = energy * (rand()-rand()) / (double)(RAND_MAX);
		V[3].id(i,j,k) = energy * (rand()-rand()) / (double)(RAND_MAX);
	}}}

	// remove mean velocities
	for (j=2; j<Ny; j++) {	// remove mean V of every XZ plane
		V[2].lyrAdd(- V[2].av(j), j);
	}
	for (i=0; i<Nx; i++) {	// remove mean U of every YZ plane
		double mean = 0.0;
		for (k=0; k<Nz; k++) mean += V[1].yMeanU(i, k) / Nz;
		for (j=1; j<Ny; j++) {
		for (k=0; k<Nz; k++) {
			V[1].id(i,j,k) -= mean;
		}}
	}
	for (k=0; k<Nz; k++) {	// remove mean W of every XY plane
		double mean = 0.0;
		for (i=0; i<Nx; i++) mean += V[3].yMeanU(i, k) / Nx;
		for (j=1; j<Ny; j++) {
		for (i=0; i<Nx; i++) {
			V[3].id(i,j,k) -= mean;
		}}
	}

	// impose parabolic profile from laminar flow (Ly is taken for 2.0)
	for (j=1; j<Ny; j++) V[1].lyrAdd(yc[j] * (2.-yc[j]), j);

	// modify flow rate
	V[1] *= 1. / V[1].bulkMeanU(); // bulk mean U = 1.0 due to non-dimensionalization
	V[3] += - V[3].bulkMeanU(); // note: boundaries are modified through using bulk functions

	// initiate boundaries with homogeneous BC
	V[1].lyrSet(V[1][1], 0).lyrMlt(- dy[0]/dy[1], 0).lyrSet(V[1][Ny-1], Ny).lyrMlt(- dy[Ny]/dy[Ny-1], Ny);
	V[3].lyrSet(V[3][1], 0).lyrMlt(- dy[0]/dy[1], 0).lyrSet(V[3][Ny-1], Ny).lyrMlt(- dy[Ny]/dy[Ny-1], Ny);
	V[2].lyrSet(0., 1).lyrSet(0., Ny);

	return *this;
}

Feld& Feld::initfrom(const Feld &feld)
/* initiate flow field (U, V, W, P, including boundaries) from existing fields */
{
	Interp(feld.V[1], V[1]).bulkInterp('U');
	Interp(feld.V[2], V[2]).bulkInterp('V');
	Interp(feld.V[3], V[3]).bulkInterp('U');
	Interp(feld.S,    S   ).bulkInterp('U');
	V[2].lyrSet(0., 0);
	return *this;
}


/***** file IO operations *****/

Feld& Feld::writeField(const char *path, int tstep, char *suffix) const
{
	char str[32];
	sprintf(str, "U%s%08i", suffix, tstep); V[1].fileIO(path, str, 'w');
	sprintf(str, "V%s%08i", suffix, tstep); V[2].fileIO(path, str, 'w');
	sprintf(str, "W%s%08i", suffix, tstep); V[3].fileIO(path, str, 'w');
	sprintf(str, "P%s%08i", suffix, tstep);    S.fileIO(path, str, 'w');
	return (Feld&)(*this);
}

Feld& Feld::readField(const char *path, int tstep, char *suffix) const
{
	char str[32];
	sprintf(str, "U%s%08i", suffix, tstep); V[1].fileIO(path, str, 'r');
	sprintf(str, "V%s%08i", suffix, tstep); V[2].fileIO(path, str, 'r');
	sprintf(str, "W%s%08i", suffix, tstep); V[3].fileIO(path, str, 'r');
	sprintf(str, "P%s%08i", suffix, tstep);    S.fileIO(path, str, 'r');
	return (Feld&)(*this);
}

void Feld::writeTecplot(const char *path, int tstep, double time) const
/* write velocity and pressure fields to ascii files readable by tecplot */
{
	FILE *fp;
	char str[1024], fn[1024];
	int i, j, k;

	sprintf(fn, "%sFIELD%08i.dat", path?path:"", tstep);

	fp = fopen(fn, "w");

	fputs("Title = \"3D instantaneous field\"\n", fp);
	fputs("variables = \"x\", \"y\", \"z\", \"u\", \"v\", \"w\", \"p\"\n", fp);
	sprintf(str, "zone t = \"%f\", i = %i, j = %i, k = %i\n", time, Nx, Nz, Ny+1);
	fputs(str, fp);

	Mesh ms(Nx,0,Nz,Lx,0,Lz);
	Scla ul(ms), vl(ms), wl(ms), pl(ms);

	for (j=0; j<=Ny; j++) {

		V.layer2CC(ul[0], vl[0], wl[0], j);
		pl = S[j];

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
	fclose(fp);

	sprintf(str, "preplot %s", fn);	system(str);
	sprintf(str, "rm %s", fn);		system(str);
}












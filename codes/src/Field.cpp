# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>

# include "Field.h"

using namespace std;



Field::Field(int *dim):
Nx(dim[0]), Ny(dim[1]), Nz(dim[2]), Nxz(dim[0]*dim[2])
{
	Nxc = (int) (Nx/2+1);
	Nxr = 2 * Nxc;
	Nxzc = Nz * Nxc;
	Nxzr = Nz * Nxr;

	// allocate memory for pointers
	u = new double [Nxz * (Ny+1)];	// Ny-1 boxes + 2 walls
	v = new double [Nxz * (Ny+1)];	// Ny layers (including 2 walls) + 1 redundant layer (index 0)
	w = new double [Nxz * (Ny+1)];	// Ny-1 boxes + 2 walls
	p = new double [Nxz * (Ny+1)];	// Ny-1 boxes + 2 walls

	uh = new double [Nxz * (Ny+1)];	// grid same as U
	vh = new double [Nxz * (Ny+1)];	// grid same as V
	wh = new double [Nxz * (Ny+1)];	// grid same as W
	dp = new double [Nxz * (Ny+1)];	// grid same as P

	ubc = new double [Nxz * 2];		// 0 -> lower wall, 1 -> upper wall
	vbc = new double [Nxz * 2];		// 0 -> lower wall, 1 -> upper wall
	wbc = new double [Nxz * 2];		// 0 -> lower wall, 1 -> upper wall

	fdp = new double [Nxzr * Ny];	// the complex party of fft should be a bit larger

	frcs = new fftw_plan [Ny-1];
	fcrs = new fftw_plan [Ny-1];

	// set up fft plans
	double *rs = (double*) &( dp[Nxz] );	// the effective domain of DP starts from the second layer, since the wall layer do not participate in computation
	double *cs = (double*) &( fdp[Nxzr] );	// the effective domain of FDP starts from the second layer, since the wall layer do not participate in computation
	
	for (int j=0; j<Ny-1; j++) {
		double *r = (double*) &( rs[Nxz * j] );
		fcmplx *c = (fcmplx*) &( cs[Nxzr * j] );
		frcs[j] = fftw_plan_dft_r2c_2d(Nz, Nx, r, c, FFTW_MEASURE);
		fcrs[j] = fftw_plan_dft_c2r_2d(Nz, Nx, c, r, FFTW_MEASURE); // note: the inverse transform (complex to real) has the side-effect of overwriting its input array
	}

}

Field::~Field()
{
	delete [] u; delete [] v; delete [] w; delete [] p;
	delete [] uh; delete [] vh; delete [] wh; delete [] dp;
	delete [] vbc; delete [] vbc; delete [] vbc;
	delete [] fdp;

	for (int j=0; j<Ny-1; j++) {
		fftw_destroy_plan(frcs[j]);
		fftw_destroy_plan(fcrs[j]);
	}
	delete [] frcs;
	delete [] fcrs;
}


void Field::initField(double energy, class Mesh *pmesh)
/* initiate flow field */
{
	int idx, i, j, k;
	double *y = pmesh->y, Ly = pmesh->Ly;

	// initiate velocities with random fluctuation
	srand(time(0));
	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		u[idx] = energy * (rand()-rand()) / (double)(RAND_MAX); // probability distribution: p(x) = ( x<0 ? x+1 : 1-x )
		v[idx] = energy * (rand()-rand()) / (double)(RAND_MAX); // 2nd order moment is 2/3, <U^2+V^2+W^2>/2 = energy
		w[idx] = energy * (rand()-rand()) / (double)(RAND_MAX); // sample space expands on the whole physical domain
	}}}

	// remove mean velocities
	for (j=2; j<Ny; j++) {	// remove mean V at every XZ plane
		this->layerAdd(v, - pmesh->layerMean(v, j), j);
	}
	for (i=0; i<Nx; i++) {	// remove mean U at every YZ plane
		double mean = 0.0;
		for (k=0; k<Nz; k++) mean += pmesh->yMeanU(u, i, k) / Nz;
		for (j=1; j<Ny; j++) {
		for (k=0; k<Nz; k++) {
			u[IDX(i,j,k)] -= mean;
		}}
	}
	for (k=0; k<Nz; k++) {	// remove mean W at every XY plane
		double mean = 0.0;
		for (i=0; i<Nx; i++) mean += pmesh->yMeanU(w, i, k) / Nx;
		for (j=1; j<Ny; j++) {
		for (i=0; i<Nx; i++) {
			w[IDX(i,j,k)] -= mean;
		}}
	}

	// impose parabolic profile from laminar flow
	for (j=1; j<Ny; j++) {
		double yh = 0.5 * (y[j] + y[j+1]);
		this->layerAdd(u, yh * (Ly-yh), j);
	}

	// modify flow rate
	this->bulkMult( u, 1.0 / pmesh->bulkMeanU(u) ); // rescale the mass flow rate to be equal to 2.0 (bulk mean U = 1.0 because of non-dimensionalization)
	this->bulkAdd( w, - pmesh->bulkMeanU(w) );
	this->bulkSet( p, 0 );
	mpg[0] = mpg[1] = mpg[2] = 0.0;

	// initiate boundary conditions
	this->bcond(0);
	// implement BC on boundaries
	this->applyBC();
}

void Field::bcond(int tstep)
/* set boundary conditions (not yet applied to velocity field) */
{
	this->layerSet(ubc, 0.0, 0);
	this->layerSet(ubc, 0.0, 1);
	this->layerSet(vbc, 0.0, 0);
	this->layerSet(vbc, 0.0, 1);
	this->layerSet(wbc, 0.0, 0);
	this->layerSet(wbc, 0.0, 1);
}

void Field::applyBC()
/* apply Dirichlet BC on velocities */
{
	this->layerCopy(u, ubc, 0 , 0);
	this->layerCopy(u, ubc, Ny, 1);
	this->layerCopy(v, vbc, 1 , 0);
	this->layerCopy(v, vbc, Ny, 1);
	this->layerCopy(w, wbc, 0 , 0);
	this->layerCopy(w, wbc, Ny, 1);
}



/***** convinent operations for whole arrays *****/

double* Field::layerCopy(double *dst, double *src, int j1, int j0)
{
	memcpy(& dst[Nxz * j1], & src[Nxz * j0], sizeof(double) * Nxz);
	return dst;
}
double* Field::layerSet(double *dst, double a, int j)
{
	int k, i, kk, jj = Nxz * j;
	for (k=0; k<Nz; k++) { kk = jj + Nx * k;
	for (i=0; i<Nx; i++) { dst[kk + i] = a;	}}
	return dst;
}
double* Field::layerAdd(double *dst, double a, int j)
{
	int k, i, kk, jj = Nxz * j;
	for (k=0; k<Nz; k++) { kk = jj + Nx * k;
	for (i=0; i<Nx; i++) { dst[kk + i] += a;	}}
	return dst;
}
double* Field::layerMult(double *dst, double a, int j)
{
	int k, i, kk, jj = Nxz * j;
	for (k=0; k<Nz; k++) { kk = jj + Nx * k;
	for (i=0; i<Nx; i++) { dst[kk + i] *= a;	}}
	return dst;
}
double* Field::layersAdd(double *dst, double *src, int j1, int j0)
{
	int k, i, kk0, kk1, jj0 = Nxz * j0, jj1 = Nxz * j1;
	for (k=0; k<Nz; k++) { kk0 = jj0 + Nx * k; kk1 = jj1 + Nx * k;
	for (i=0; i<Nx; i++) { dst[kk1 + i] += src[kk0 + i];	}}
	return dst;
}
double* Field::layersMult(double *dst, double *src, int j1, int j0)
{
	int k, i, kk0, kk1, jj0 = Nxz * j0, jj1 = Nxz * j1;
	double *src1 = src ? src : dst; // !!! CAUTION: when src1 == dst, all write to dst memory also changes src memory
	for (k=0; k<Nz; k++) { kk0 = jj0 + Nx * k; kk1 = jj1 + Nx * k;
	for (i=0; i<Nx; i++) { dst[kk1 + i] *= src1[kk0 + i];	}}
	return dst;
}

double* Field::bulkCopy(double *dst, double *src)
{
	memcpy(& dst[Nxz], & src[Nxz], sizeof(double) * Nxz * (Ny-1));
	return dst;
}
double* Field::bulkSet(double *dst, double a)
{
	for (int j=1; j<Ny; j++)	this->layerSet(dst, a, j);
	return dst;
}
double* Field::bulkAdd(double *dst, double a)
{
	for (int j=1; j<Ny; j++)	this->layerAdd(dst, a, j);
	return dst;
}
double* Field::bulkMult(double *dst, double a)
{
	for (int j=1; j<Ny; j++)	this->layerMult(dst, a, j);
	return dst;
}


/***** file input & output operations *****/

void Field::writeField(char *path, int tstep)
{
	double *ptrs[4] = {u, v, w, p};
	char names[4][32] = {"U", "V", "W", "P"};
	for (int n=0; n<4; n++) {
		sprintf(names[n], "%s%08i", names[n], tstep);
		this->fieldIO(path, ptrs[n], names[n], 'w');
	}
}

void Field::readField(char *path, int tstep)
{
	double *ptrs[4] = {u, v, w, p};
	char names[4][32] = {"U", "V", "W", "P"};
	for (int n=0; n<4; n++) {
		sprintf(names[n], "%s%08i", names[n], tstep);
		this->fieldIO(path, ptrs[n], names[n], 'r');
	}
}

void Field::fieldIO(char *path, double *ptr, char *name, char mode)
/* read & write velocities and pressure from & to binary files */
{
	FILE *fp;
	char str[1024];

	sprintf(str, "%s%s.bin", path?path:"", name);

	fp = fopen( str, (mode == 'w' ? "wb" : "rb") );

		// write domain information at the beginning
		if (mode == 'w') {
			fwrite(& Nx, sizeof(int), 1, fp);
			fwrite(& Ny, sizeof(int), 1, fp);
			fwrite(& Nz, sizeof(int), 1, fp);
		}
		// data begin after the info section
		fseek(fp, sizeof(double) * Nxz, SEEK_SET);
		if (mode == 'w') fwrite(ptr, sizeof(double) * Nxz, Ny+1, fp); // note: first layer (j=0) of V is redundant
		if (mode == 'r') fread (ptr, sizeof(double) * Nxz, Ny+1, fp);

	fclose(fp);
}

void Field::writeTecplot(char *path, int tstep, double time, class Mesh *pmesh)
/* write velocity and pressure fields to ascii files readable by tecplot */
{
	FILE *fp;
	char str[1024];
	int i, j, k, idx, ip, jp, kp;
	double v1, v2;

	sprintf(str, "%sFIELD%08i.plt", path?path:"", tstep);

	fp = fopen(str, "w");

		fputs("Title = \"3D instantaneous field\"\n", fp);
		fputs("variables = \"x\", \"z\", \"y\", \"u\", \"v\", \"w\", \"p\"\n", fp);
		sprintf(str, "zone t = \"%f\", i = %i, j = %i, k = %i\n", time, Ny+1, Nz, Nx);
		fputs(str, fp);

		for (i=0; i<Nx; i++) {
		for (k=0; k<Nz; k++) {
		for (j=0; j<=Ny;j++) {
			idx = IDX(i,j,k);
			ip = IDX(pmesh->ipa[i],j,k);
			jp = IDX(i,j+1,k);
			kp = IDX(i,j,pmesh->kpa[k]);

			v1 = j==0 ? v[jp] : v[idx];
			v2 = j==Ny ? v[idx] : v[jp];
			// all interpolate to cell centers
			sprintf(str, "%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
				pmesh->dx * (i + 0.5),
				pmesh->dz * (k + 0.5),
				pmesh->yc[j],
				0.5 * ( u[idx] + u[ip] ),
				0.5 * ( v1 + v2 ),
				0.5 * ( w[idx] + w[kp] ),
				p[idx]	);
			fputs(str, fp);
		}}}

	fclose(fp);
}




/***** debug functions *****/

void Field::debug_AsciiOutput(int tstep, double *ptr, int j1, int j2, string filepre)
/* write the fields in ascii files for check */
{
	int i, j, k;
	char str[1024];	sprintf(str, "debug/%s%08i.txt", filepre.c_str(), tstep);
	FILE *fp = fopen(str, "w");

	for (j=j1;j<=j2;j++) {	fputc('\n', fp);
	for (k=0; k<Nz; k++) {	fputc('\n', fp);
	for (i=0; i<Nx; i++) {
		sprintf(str, "%.6f\t", ptr[IDX(i,j,k)]);
		fputs(str, fp);
	}}}
	fclose(fp);
}

void Field::debug_AsciiOutputF(int tstep, int j1, int j2)
/* write the Fourier transformed fields in ascii files for check */
{
	int i, j, k, idx, nxc = (int) (Nx/2+1), nxr = 2 * nxc;
	char str[1024];
	FILE *fp1, *fp2;

	sprintf(str, "debug/FR%08i.txt", tstep);	fp1 = fopen(str, "w");
	sprintf(str, "debug/FI%08i.txt", tstep);	fp2 = fopen(str, "w");

	for (j=j1;j<=j2;j++) {	fputc('\n', fp1);	fputc('\n', fp2);
	for (k=0; k<Nz; k++) {	fputc('\n', fp1);	fputc('\n', fp2);
	for (i=0; i<nxc;i++) {
		idx = Nz*nxr * j + nxr * k + 2*i;
		sprintf(str, "%.6f\t", fdp[idx]);	fputs(str, fp1);
		sprintf(str, "%.6f\t", fdp[idx+1]);	fputs(str, fp2);
	}}}
	fclose(fp1);
	fclose(fp2);
}

void Field::debug_Output(int tstep)
{
	this->debug_AsciiOutput(tstep, u, 0, Ny, "U");
	this->debug_AsciiOutput(tstep, v, 1, Ny, "V");
	this->debug_AsciiOutput(tstep, w, 0, Ny, "W");
	this->debug_AsciiOutput(tstep, p, 1, Ny-1, "P");

	this->debug_AsciiOutput(tstep, uh, 0, Ny, "Uh");
	this->debug_AsciiOutput(tstep, vh, 1, Ny, "Vh");
	this->debug_AsciiOutput(tstep, wh, 0, Ny, "Wh");
	this->debug_AsciiOutput(tstep, dp, 1, Ny-1, "DP");

	this->debug_AsciiOutput(tstep, ubc, 0, 1, "Ubc");
	this->debug_AsciiOutput(tstep, vbc, 0, 1, "Vbc");
	this->debug_AsciiOutput(tstep, wbc, 0, 1, "Wbc");

	this->debug_AsciiOutputF(tstep, 1, Ny-1);
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








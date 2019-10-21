# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <cmath>

# include "Field.h"

using namespace std;



Field::Field(int dim[3]):
Nx(dim[0]), Ny(dim[1]), Nz(dim[2]), Nxz(dim[0]*dim[2]),
Nxc(Nx/2+1), Nxr(2*Nxc), Nxzc(Nz*Nxc), Nxzr(Nz*Nxr)
{
	// Nxc = (int) (Nx/2+1);
	// Nxr = 2 * Nxc;
	// Nxzc = Nz * Nxc;
	// Nxzr = Nz * Nxr;

	// allocate memory for pointers
	u = new double [Nxz * (Ny+1)];	// Ny-1 boxes + 2 walls
	v = new double [Nxz * (Ny+1)];	// Ny layers (including 2 walls) + 1 redundant layer (index 0)
	w = new double [Nxz * (Ny+1)];	// Ny-1 boxes + 2 walls
	p = new double [Nxz * (Ny+1)];	// Ny-1 boxes + 2 walls

	uh = new double [Nxz * (Ny+1)];	// grid same as U
	vh = new double [Nxz * (Ny+1)];	// grid same as V
	wh = new double [Nxz * (Ny+1)];	// grid same as W
	dp = new double [Nxz * (Ny+1)];	// grid same as P
	nu = new double [Nxz * (Ny+1)];	// grid same as P

	ubc = new double [Nxz * 2];		// 0 -> lower wall, 1 -> upper wall
	vbc = new double [Nxz * 2];		// 0 -> lower wall, 1 -> upper wall
	wbc = new double [Nxz * 2];		// 0 -> lower wall, 1 -> upper wall

	fdp = new double [Nxzr * (Ny+1)];	// the complex party of fft should be a bit larger

	U[0] = UP[0] = u; UH[0] = UPH[0] = uh; UBC[0] = ubc;
	U[1] = UP[1] = v; UH[1] = UPH[1] = vh; UBC[1] = vbc;
	U[2] = UP[2] = w; UH[2] = UPH[2] = wh; UBC[2] = wbc;
	P[0] = UP[3] = p; DP[0] = UPH[3] = dp;
	P[1] = mpg; DP[1] = fdp;

	frcs = new fftw_plan [Ny+1];
	fcrs = new fftw_plan [Ny+1];

	for (int j=0; j<=Ny; j++) {	// note: the effective domain of DP actually starts from the second layer, since the wall layers do not participate in computation
		double *r = (double*) &( dp[Nxz * j] );
		fcmplx *c = (fcmplx*) &( fdp[Nxzr * j] );
		frcs[j] = fftw_plan_dft_r2c_2d(Nz, Nx, r, c, FFTW_MEASURE);
		fcrs[j] = fftw_plan_dft_c2r_2d(Nz, Nx, c, r, FFTW_MEASURE); // note: the inverse transform (complex to real) has the side-effect of overwriting its input array
	}

}

Field::~Field()
{
	delete [] u; delete [] v; delete [] w; delete [] p;
	delete [] uh; delete [] vh; delete [] wh; delete [] dp;
	delete [] ubc; delete [] vbc; delete [] wbc;
	delete [] fdp; delete [] nu;
	for (int j=0; j<=Ny; j++) {
		fftw_destroy_plan(frcs[j]);
		fftw_destroy_plan(fcrs[j]);
	}
	delete [] frcs;
	delete [] fcrs;
}


void Field::initField(class Field *pf0, class Mesh *pm0, class Mesh *pm)
/* initiate flow field from given fields, interpolated to the current grid */
{
	int i, j, k, idx, j0, idx0, jm0;
	int *i0 = new int [Nx], *k0 = new int [Nz];
	int Nx0 = pm0->Nx, Ny0 = pm0->Ny, Nz0 = pm0->Nz;
	double *y = pm->y, *yc = pm->yc, *y0 = pm0->y, *yc0 = pm0->yc;

	// nearest interpolation for x and z directions in Fourier space
	// find out nearest point for every wavenumber k_z
	for (k=0; k<Nz; k++) {	k0[k] = 0;
	for (int k1=1; k1<Nz0; k1++) {
		if (	fabs( pm->kz(k) - pm0->kz(k1)	)
			<	fabs( pm->kz(k) - pm0->kz(k0[k])	)	)	k0[k] = k1;
	}}
	// find out nearesr point for every wavenumber k_x
	for (i=0; i<Nx; i++) {	i0[i] = 0;
	for (int i1=1; i1<Nx0; i1++) {
		if (	fabs( pm->kx(i) - pm0->kx(i1)	)
			<	fabs( pm->kx(i) - pm0->kx(i0[i])	)	)	i0[i] = i1;
	}}

	for (int n=0; n<4; n++) {

		pf0->bulkCopy(pf0->dp, pf0->UP[n]);
		pf0->fft();

		// linear interpolations for y direction
		for (j=0; j<=Ny; j++) {
			// find out the range (j0-1 ~ j0) where j fall
			for (j0=0; j0<=Ny0; j0++) { if (yc0[j0] > yc[j]) break; }
			if (n == 1) { j0 = j==0 ? 0 : (j==1 ? 2 : j0); }

			for (k=0; k<Nz; k++) {
			for (i=0; i<Nxc; i++) {

				idx = this->IDXF(i,j,k);
				idx0 = pf0->IDXF(i0[i],j0,k0[k]);
				jm0 = pf0->IDXF(i0[i],j0-1,k0[k]);

				// values at both ends are extended for out-of-range points
				if (j0 == 0) {
					fdp[idx] = pf0->fdp[idx0];		// real part
					fdp[idx+1] = pf0->fdp[idx0+1];	// imaginary part
				}
				else if (j0 == Ny0+1) {
					fdp[idx] = pf0->fdp[jm0];
					fdp[idx+1] = pf0->fdp[jm0+1];
				}
				// linear interpolation for in-range points
				else {
					fdp[idx] =
						(yc0[j0] - yc[j]) / ( yc0[j0] - yc0[j0-1] ) * pf0->fdp[jm0]
					+	(yc[j]-yc0[j0-1]) / ( yc0[j0] - yc0[j0-1] ) * pf0->fdp[idx0];

					fdp[idx+1] =
						(yc0[j0] - yc[j]) / ( yc0[j0] - yc0[j0-1] ) * pf0->fdp[jm0+1]
					+	(yc[j]-yc0[j0-1]) / ( yc0[j0] - yc0[j0-1] ) * pf0->fdp[idx0+1];
				}
			}
		}}
		this->ifft();
		this->bulkCopy( this->UP[n], this->bulkMult(this->dp, 1.0/pm0->Nxz*this->Nxz) );
	}

	// 3 mean pressure gradients are initiated to 0, and will be compensated by dmpg next step
	mpg[0] = mpg[1] = mpg[2] = 0.0;

	// // implement BC on boundaries
	// this->applyBC();

	printf("\nFlow fields initiated from existing fields.\n");

	delete [] i0; delete [] k0;
}

void Field::initField(double energy, class Mesh *pmesh)
/* initiate flow field from laminar with random fluctions */
{
	int idx, i, j, k;
	double *y = pmesh->y, *yc = pmesh->yc, Ly = pmesh->Ly;

	// initiate velocities with random fluctuations
	srand(time(0));
	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		u[idx] = energy * (rand()-rand()) / (double)(RAND_MAX); // probability distribution: p(x) = ( x<0 ? x+1 : 1-x )
		if (j>1) v[idx] = energy * (rand()-rand()) / (double)(RAND_MAX); // 2nd order moment is 2/3, <U^2+V^2+W^2>/2 = energy
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
	for (j=1; j<Ny; j++)	this->layerAdd(u, yc[j] * (Ly-yc[j]), j);

	// modify flow rate // note: boundary will be modified through using bulk functions
	this->bulkMult( u, 1.0 / pmesh->bulkMeanU(u) ); // rescale the mass flow rate to be equal to 2.0 (bulk mean U = 1.0 because of non-dimensionalization)
	this->bulkAdd( w, - pmesh->bulkMeanU(w) );

	// set initial pressure and its mean gradients in 3 directions
	this->bulkSet( p, 0 );
	mpg[0] = mpg[1] = mpg[2] = 0.0;

	// implement BC on boundaries
	this->bcond(0);
	this->applyBC();

	printf("\nFlow fields initiated from laminar.\n");
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

void Field::applyBC(double dt)
/* apply Dirichlet BC on velocities */
{
	double *ul = new double [Nxz];
	double *vl = new double [Nxz];
	double *wl = new double [Nxz];
	
	this->layerMult( layerCopy(ul,u,0,0), -1 );
	this->layerMult( layerCopy(vl,v,0,1), -1 );
	this->layerMult( layerCopy(wl,w,0,0), -1 );
	this->layerMult( layersAdd(ul,ubc,0,0), 1.0 / dt );
	this->layerMult( layersAdd(vl,vbc,0,0), 1.0 / dt );
	this->layerMult( layersAdd(wl,wbc,0,0), 1.0 / dt );
	this->layerCopy(uh, ul, 0, 0);
	this->layerCopy(vh, vl, 1, 0);
	this->layerCopy(wh, wl, 0, 0);
	
	this->layerMult( layerCopy(ul,u,0,Ny), -1 );
	this->layerMult( layerCopy(vl,v,0,Ny), -1 );
	this->layerMult( layerCopy(wl,w,0,Ny), -1 );
	this->layerMult( layersAdd(ul,ubc,0,1), 1.0 / dt );
	this->layerMult( layersAdd(vl,vbc,0,1), 1.0 / dt );
	this->layerMult( layersAdd(wl,wbc,0,1), 1.0 / dt );
	this->layerCopy(uh, ul, Ny, 0);
	this->layerCopy(vh, vl, Ny, 0);
	this->layerCopy(wh, wl, Ny, 0);

	this->layerCopy(u, ubc, 0, 0);
	this->layerCopy(v, vbc, 1, 0);
	this->layerCopy(w, wbc, 0, 0);
	this->layerCopy(u, ubc, Ny, 1);
	this->layerCopy(v, vbc, Ny, 1);
	this->layerCopy(w, wbc, Ny, 1);

	delete [] ul; delete [] vl; delete [] wl;
}




// void Field::bodyForce(int bftype)
// {
// 	int i, j, k, idx;

// 	switch (bftype) {
// 		case 1:
// 			double um, *usm = new double [Nx];
// 			double vm, *vsm = new double [Nx];
// 			double uhm, *uhsm = new double [Nx];
// 			double vhm, *vhsm = new double [Nx];

// 			for (j=1; j<Ny; j++) {

// 				for (i=0; i<Nx; i++) {
// 					usm[i] = vsm[i] = 0.0;
// 					uhsm[i] = vhsm[i] = 0.0;
// 					for (k=0; k<Nz; k++) {
// 						idx = IDX(i,j,k);
// 						usm[i] += u[idx] / Nz;
// 						vsm[i] += v[idx] / Nz;
// 						uhsm[i] += uh[idx] / Nz;
// 						vhsm[i] += vh[idx] / Nz;
// 				}}

// 				um = vm = 0.0;
// 				uhm = vhm = 0.0;
// 				for (i=0; i<Nx; i++) {
// 					um += usm[i] / Nx;
// 					vm += vsm[i] / Nx;
// 					uhm += uhsm[i] / Nx;
// 					vhm += vhsm[i] / Nx;
// 				}

// 				for (k=0; k<Nz; k++) {
// 				for (i=0; i<Nx; i++) {
// 					idx = IDX(i,j,k);
// 					u[idx] -= usm[i] - um;
// 					if (j>1) v[idx] -= vsm[i] - vm;
// 					uh[idx] -= uhsm[i] - uhm;
// 					if (j>1) vh[idx] -= vhsm[i] - vhm;
// 				}}
// 			}
// 		break;

// 		case 2:
// 		break;
// 	}
// }


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
{	// note: bulk functions will change the boundary
	memcpy(dst, src, sizeof(double) * Nxz * (Ny+1));
	return dst;
}
double* Field::bulkSet(double *dst, double a)
{	// note: bulk functions will change the boundary
	for (int j=0; j<=Ny; j++)	this->layerSet(dst, a, j);
	return dst;
}
double* Field::bulkAdd(double *dst, double a)
{	// note: bulk functions will change the boundary
	for (int j=0; j<=Ny; j++)	this->layerAdd(dst, a, j);
	return dst;
}
double* Field::bulkMult(double *dst, double a)
{	// note: bulk functions will change the boundary
	for (int j=0; j<=Ny; j++)	this->layerMult(dst, a, j);
	return dst;
}

double* Field::removeSpanMean(double *dst, int j)
{
	int i, k;
	double qm = 0.0, *qsm = new double [Nx];

	for (i=0; i<Nx; i++) {	qsm[i] = 0.0;
	for (k=0; k<Nz; k++) {	qsm[i] += dst[IDX(i,j,k)] / Nz;	}}

	for (i=0; i<Nx; i++)	qm += qsm[i] / Nx;

	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {	dst[IDX(i,j,k)] -= qsm[i] - qm;	}}

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

void Field::writeFieldDt(char *path, int tstep)
{
	double *ptrs[4] = {uh, vh, wh, dp};
	char names[4][32] = {"UT", "VT", "WT", "PT"};
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
	char str[1024], fn[1024];
	int i, j, k, idx, ip, jp, kp;
	double v1, v2;

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

	sprintf(str, "preplot %s", fn);	system(str);
	sprintf(str, "rm %s", fn);		system(str);
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








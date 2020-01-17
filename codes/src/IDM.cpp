# include <iostream>
# include <math.h>
# include <omp.h>

# include "IDM.h"
# include "Matrix.h"

using namespace std;


/***** configuration *****/

int IDM::ompset(int n)
{
	int nprocs = omp_get_num_procs(), ngrids, nthrds;
	if (n > 0) nthrds = (n > nprocs ? nprocs : n);
	else {
	/* automatically decide number of OpenMP threads based on grid number */
		ngrids = log2(Nx*Ny*Nz) - 16;
		ngrids = (ngrids < 1 ? 1 : ngrids);
		nprocs = (nprocs > 32 ? 32 : nprocs);
		nthrds = (ngrids > nprocs ? nprocs : ngrids);
	}
	omp_set_num_threads(nthrds);
	return nthrds;
}


/***** subroutines for computation *****/

void IDM::update(Feld &FLD, Feld &FLDH, double dt)
/* project from U^* to U^n+1 using DP */
{
	int i, j, k, idx, im, jm, km;

	double *u,*v,*w,*p;      FLD.ptrGet(u,v,w,p);
	double *uh,*vh,*wh,*dp; FLDH.ptrGet(uh,vh,wh,dp);
	Scla &DP = FLDH.S;

	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

		// store time derivative in intermediate variables
		uh[idx] = (uh[idx] - u[idx]) / dt - (dp[idx] - dp[im]) / dx; if (j>1) // act on next line
		vh[idx] = (vh[idx] - v[idx]) / dt - (dp[idx] - dp[jm]) / h[j];
		wh[idx] = (wh[idx] - w[idx]) / dt - (dp[idx] - dp[km]) / dz;
		// update fields
		u[idx] += uh[idx] * dt; if (j>1) // act on next line
		v[idx] += vh[idx] * dt;
		w[idx] += wh[idx] * dt;
		p[idx] += dp[idx];
	}}}

	for (j=1; j<Ny; j++) DP.lyrMlt(1./dt, j);
}


void IDM::meanpg(Vctr &U, double mpg[3], double dt)
/* solve the increment of mean pressure gradient at n+1/2 step, given mass flow rate at n+1 step */
{
	// solve mean pressure gradient increment given streamwise flow rate 2.0 and spanwise flow rate 0
	double dmpg1 = (U[1].bulkMeanU() - 1.) / dt;
	double dmpg3 =  U[3].bulkMeanU()       / dt;
	// update the mean pressure gradient
	mpg[0] += dmpg1; // here mpg is passed in by a copy of its initial address
	mpg[2] += dmpg3; // changes made to elements will act on the real array
	// complement the mean pressure gradient increment that was not included in the velocity update step
	for (int j=1; j<Ny; j++) {
		U[1].lyrAdd(- dt * dmpg1, j);
		U[3].lyrAdd(- dt * dmpg3, j);
	}
}
 

void IDM::pressBD(Scla &P, Scla &DP, const Scla &PBC)
/* boundary pressure never participate in computation, set to 0 */
{
	 P.lyrSet(0., 0).lyrSet(0., Ny);
	DP.lyrSet(0., 0).lyrSet(0., Ny);
}

void IDM::veldtBD(Vctr &UH, const Vctr &U, const Vctr &UBC, double dt)
/* modify boundary of velocity-time-derivative with given BC */
{
	Mesh ms(Nx,0,Nz,Lx,0,Lz); Scla ql(ms);

	UH[1].lyrSet( (( (ql=UBC[1][0]) -= U[1][0] ) *= 1./dt)[0], 0 );
	UH[2].lyrSet( (( (ql=UBC[2][0]) -= U[2][1] ) *= 1./dt)[0], 1 );
	UH[3].lyrSet( (( (ql=UBC[3][0]) -= U[3][0] ) *= 1./dt)[0], 0 );

	UH[1].lyrSet( (( (ql=UBC[1][1]) -= U[1][Ny] ) *= 1./dt)[0], Ny );
	UH[2].lyrSet( (( (ql=UBC[2][1]) -= U[2][Ny] ) *= 1./dt)[0], Ny );
	UH[3].lyrSet( (( (ql=UBC[3][1]) -= U[3][Ny] ) *= 1./dt)[0], Ny );

	ms.freeall();
}

void IDM::velocBD(Vctr &U, const Vctr &UBC)
/* apply Dirichlet BC on velocities */
{
	U[1].lyrSet(UBC[1][0], 0).lyrSet(UBC[1][1], Ny);
	U[2].lyrSet(UBC[2][0], 1).lyrSet(UBC[2][1], Ny);
	U[3].lyrSet(UBC[3][0], 0).lyrSet(UBC[3][1], Ny);

	// extrapolate U to virtual boundary using UBC (real boundary) at new time step
	// for (int k=0; k<Nz; k++) {
	// for (int i=0; i<Nx; i++) {
	// 	U[1].id(i,0,k) = 2.*h[1]/dy[1] * UBC[1].id(i,0,k) - dy[0]/dy[1] * U[1].id(i,1,k);
	// 	U[3].id(i,0,k) = 2.*h[1]/dy[1] * UBC[3].id(i,0,k) - dy[0]/dy[1] * U[3].id(i,1,k);
	// 	U[2].id(i,1,k) = UBC[2].id(i,0,k);

	// 	U[1].id(i,Ny,k) = 2.*h[Ny]/dy[Ny-1] * UBC[1].id(i,1,k) - dy[Ny]/dy[Ny-1] * U[1].id(i,Ny-1,k);
	// 	U[3].id(i,Ny,k) = 2.*h[Ny]/dy[Ny-1] * UBC[3].id(i,1,k) - dy[Ny]/dy[Ny-1] * U[3].id(i,Ny-1,k);
	// 	U[2].id(i,Ny,k) = UBC[2].id(i,1,k);
	// }}
}



/***** intermediate velocity computation *****/

void IDM::urhs1(double *ruh, const Feld &FLD, const Feld &VIS, double *fbx)
/* compute right hand side R_1 for intermediate velocity at all non-wall grid points */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, imjp, imkp, imjm, imkm, jup, jum;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double u1, u2, v1, v2, w1, w2;
	double api, aci, ami, apj, acj, amj, apk, ack, amk, mbcu, mbcd;
	double l11un, l12vn, l13wn, m11un, m12vn, m13wn, pressg, mbc=0;

	double *u,*v,*w,*p;        FLD.ptrGet(u,v,w,p);
	double *nux,*nuy,*nuz,*nu; VIS.ptrGet(nux,nuy,nuz,nu);

	# pragma omp for
	for (j=1; j<Ny; j++) { jup = (j!=Ny-1); jum = (j!=1); // indicate the secondary boundary
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		imjp = IDX(ima[i],j+1,k); imkp = IDX(ima[i],j,kpa[k]);
		imjm = IDX(ima[i],j-1,k); imkm = IDX(ima[i],j,kma[k]);
		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

		// interpolate viscosity and velocity at step n to the position needed
		vis1 = nu[im];   vis2 = nu[idx];
		vis3 = nuz[idx]; vis4 = nuz[jp];
		vis5 = nuy[idx]; vis6 = nuy[kp];

		u2 = 0.5 * (u[idx]+ u[ip]);
		u1 = 0.5 * (u[idx]+ u[im]);
		v2 = 0.5 * (v[jp] + v[imjp]);
		v1 = 0.5 * (v[idx]+ v[im]);
		w2 = 0.5 * (w[kp] + w[imkp]);
		w1 = 0.5 * (w[idx]+ w[im]);

		// viscous terms
		api = 1.0/dx2 * vis2;
		aci =-1.0/dx2 *(vis2+vis1);
		ami = 1.0/dx2 * vis1;
		apj = 0.5/dy[j] * vis4/h[j+1];
		acj =-0.5/dy[j] *(vis4/h[j+1]+vis3/h[j]);
		amj = 0.5/dy[j] * vis3/h[j];
		apk = 0.5/dz2 * vis6;
		ack =-0.5/dz2 *(vis6+vis5);
		amk = 0.5/dz2 * vis5;

		l11un =	api*u[ip] + aci*u[idx] + ami*u[im]
			+	apj*u[jp] + acj*u[idx] + amj*u[jm]
			+	apk*u[kp] + ack*u[idx] + amk*u[km];
		l12vn = ( vis4 * (v[jp]-v[imjp]) - vis3 * (v[idx]-v[im]) ) / (2.0*dx*dy[j]);
		l13wn = ( vis6 * (w[kp]-w[imkp]) - vis5 * (w[idx]-w[im]) ) / (2.0*dx*dz);

		// non-linear terms
		// m11un
		api = 0.5/dx * u2     - api;
		aci = 0.5/dx *(u2-u1) - aci;
		ami =-0.5/dx * u1     - ami;
		apj = 0.25/h[j+1]* v2                                   - apj;
		acj = 0.25/dy[j] *(v2*dy[j+1]/h[j+1] - v1*dy[j-1]/h[j]) - acj;
		amj =-0.25/h[j]  * v1                                   - amj;
		apk = 0.25/dz * w2     - apk;
		ack = 0.25/dz *(w2-w1) - ack;
		amk =-0.25/dz * w1     - amk;

		// mbcd = amj * ubc[IDX(i,0,k)];
		// mbcu = apj * ubc[IDX(i,1,k)];
		// mbc = (jum-1) * mbcd + (jup-1) * mbcu;
		apj *= jup;
		amj *= jum;
		m11un =	api*u[ip] + aci*u[idx] + ami*u[im]
			+	apj*u[jp] + acj*u[idx] + amj*u[jm]
			+	apk*u[kp] + ack*u[idx] + amk*u[km];

		// m12vn
		u2 = ( u[idx]*dy[j+1] + u[jp]*dy[j] ) / (2.0*h[j+1]);
		u1 = ( u[idx]*dy[j-1] + u[jm]*dy[j] ) / (2.0*h[j]);

		// mbcd = .5/dy[1]    * (-u1 * (vbc[IDX(i,0,k)]+vbc[IDX(ima[i],0,k)])/2. + vis3 * (vbc[IDX(i,0,k)]-vbc[IDX(ima[i],0,k)])/dx );
		// mbcu = .5/dy[Ny-1] * ( u2 * (vbc[IDX(i,1,k)]+vbc[IDX(ima[i],1,k)])/2. - vis4 * (vbc[IDX(i,1,k)]-vbc[IDX(ima[i],1,k)])/dx );
		// mbc += (jum-1) * mbcd + (jup-1) * mbcu;
		u2 *= jup;	vis4 *= jup;
		u1 *= jum;	vis3 *= jum;
		m12vn = .5/dy[j] * ( u2*v2 - u1*v1 - ( vis4 * (v[jp]-v[imjp]) - vis3 * (v[idx]-v[im]) )/dx );

		// m13wn
		u2 = 0.5 * ( u[idx] + u[kp] );
		u1 = 0.5 * ( u[idx] + u[km] );
		m13wn = (u2*w2 - u1*w1) / (2.0*dz) - l13wn;

		// pressure gradient term
		pressg = (p[idx] - p[im]) / dx;

		// vis3 = nuz[idx];
		// vis4 = nuz[jp];
		// v2 = 0.5 * (v[jp] + v[imjp]);
		// v1 = 0.5 * (v[idx]+ v[im]);
		// amj = -0.25/h[j]  * v1 - 0.5/dy[j] * vis3/h[j];
		// apj = 0.25/h[j+1]* v2 - 0.5/dy[j] * vis4/h[j+1];
		// u2 = ( u[idx]*dy[j+1] + u[jp]*dy[j] ) / (2.0*h[j+1]);
		// u1 = ( u[idx]*dy[j-1] + u[jm]*dy[j] ) / (2.0*h[j]);
		// mbcd = amj * ubc[IDX(i,0,k)] + .5/dy[1]    * (-u1 * (vbc[IDX(i,0,k)]+vbc[IDX(ima[i],0,k)])/2. + vis3 * (vbc[IDX(i,0,k)]-vbc[IDX(ima[i],0,k)])/dx );
		// mbcu = apj * ubc[IDX(i,1,k)] + .5/dy[Ny-1] * ( u2 * (vbc[IDX(i,1,k)]+vbc[IDX(ima[i],1,k)])/2. - vis4 * (vbc[IDX(i,1,k)]-vbc[IDX(ima[i],1,k)])/dx );
		// mbc = (jum-1) * mbcd + (jup-1) * mbcu;

		// R_1 without boundary modification
		ruh[idx] = (l11un + l12vn + l13wn)
		         - (m11un + m12vn + m13wn)
		         - pressg + fbx[idx] + mbc;
	}}}
}

void IDM::urhs2(double *rvh, const Feld &FLD, const Feld &VIS, double *fby)
/* compute right hand side R_2 for intermediate velocity at all non-wall grid points */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, ipjm, jmkp, imjm, jmkm, jup, jum;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double u1, u2, v1, v2, w1, w2;
	double api, aci, ami, apj, acj, amj, apk, ack, amk, mbcd, mbcu;
	double l21un, l22vn, l23wn, m21un, m22vn, m23wn, pressg, mbc=0;

	double *u,*v,*w,*p;        FLD.ptrGet(u,v,w,p);
	double *nux,*nuy,*nuz,*nu; VIS.ptrGet(nux,nuy,nuz,nu);

	# pragma omp for
	for (j=2; j<Ny; j++) { jup = (j!=Ny-1); jum = (j!=2);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		ipjm = IDX(ipa[i],j-1,k); jmkp = IDX(i,j-1,kpa[k]);
		imjm = IDX(ima[i],j-1,k); jmkm = IDX(i,j-1,kma[k]);
		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

		// interpolate viscosity and velocity at step n to the position needed
		vis1 = nuz[idx]; vis2 = nuz[ip];
		vis3 = nu[jm];   vis4 = nu[idx];
		vis5 = nux[idx]; vis6 = nux[kp];

		u2 = ( u[ip]*dy[j-1] + u[ipjm]*dy[j] ) / (2.0*h[j]);
		u1 = ( u[idx]*dy[j-1] + u[jm]*dy[j] ) / (2.0*h[j]);
		v2 = 0.5 * ( v[idx] + v[jp] );
		v1 = 0.5 * ( v[idx] + v[jm] );
		w2 = ( w[kp]*dy[j-1] + w[jmkp]*dy[j] ) / (2.0*h[j]);
		w1 = ( w[idx]*dy[j-1] + w[jm]*dy[j] ) / (2.0*h[j]);

		// viscous terms
		api = 0.5/dx2 * vis2;
		aci =-0.5/dx2 *(vis2+vis1);
		ami = 0.5/dx2 * vis1;
		apj = 1.0/h[j] * vis4/dy[j];
		acj =-1.0/h[j] *(vis4/dy[j]+vis3/dy[j-1]);
		amj = 1.0/h[j] * vis3/dy[j-1];
		apk = 0.5/dz2 * vis6;
		ack =-0.5/dz2 *(vis6+vis5);
		amk = 0.5/dz2 * vis5;

		l21un = ( vis2 * (u[ip]-u[ipjm]) - vis1 * (u[idx]-u[jm]) ) / (2.0*h[j]*dx);
		l22vn =	api*v[ip] + aci*v[idx] + ami*v[im]
			+	apj*v[jp] + acj*v[idx] + amj*v[jm]
			+	apk*v[kp] + ack*v[idx] + amk*v[km];
		l23wn = ( vis6 * (w[kp]-w[jmkp]) - vis5 * (w[idx]-w[jm]) ) / (2.0*h[j]*dz);

		// non-linear terms
		// m22vn
		api = 0.25/dx * u2     - api;
		aci = 0.25/dx *(u2-u1) - aci;
		ami =-0.25/dx * u1     - ami;
		apj = 0.5/h[j]* v2     - apj;
		acj = 0.5/h[j]*(v2-v1) - acj;
		amj =-0.5/h[j]* v1     - amj;
		apk = 0.25/dz * w2     - apk;
		ack = 0.25/dz *(w2-w1) - ack;
		amk =-0.25/dz * w1     - amk;

		// mbcd = amj * vbc[IDX(i,0,k)];
		// mbcu = apj * vbc[IDX(i,1,k)];
		// mbc = (jum-1) * mbcd + (jup-1) * mbcu;
		apj *= jup;
		amj *= jum;
		m22vn =	api*v[ip] + aci*v[idx] + ami*v[im]
			+	apj*v[jp] + acj*v[idx] + amj*v[jm]
			+	apk*v[kp] + ack*v[idx] + amk*v[km];

		// m21un
		v2 = 0.5 * ( v[idx] + v[ip] );
		v1 = 0.5 * ( v[idx] + v[im] );
		m21un = (v2*u2 - v1*u1) / (2.0*dx) - l21un;

		// m23wn
		v2 = 0.5 * ( v[idx] + v[kp] );
		v1 = 0.5 * ( v[idx] + v[km] );
		m23wn = (v2*w2 - v1*w1) / (2.0*dz) - l23wn;

		// pressure gradient term
		pressg = (p[idx] - p[jm]) / h[j];

		// v2 = 0.5 * ( v[idx] + v[jp] );
		// v1 = 0.5 * ( v[idx] + v[jm] );
		// amj = -0.5/h[j]* v1 - 1.0/h[j] * vis3/dy[j-1];
		// apj = 0.5/h[j]* v2 - 1.0/h[j] * vis4/dy[j];
		// mbcd = amj * vbc[IDX(i,0,k)];
		// mbcu = apj * vbc[IDX(i,1,k)];
		// mbc = (jum-1) * mbcd + (jup-1) * mbcu;

		// R_2 without boundary modification
		rvh[idx] = (l21un + l22vn + l23wn)
		         - (m21un + m22vn + m23wn)
		         - pressg + fby[idx] + mbc;
	}}}
}

void IDM::urhs3(double *rwh, const Feld &FLD, const Feld &VIS, double *fbz)
/* compute right hand side R_3 for intermediate velocity at all non-wall grid points */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, ipkm, jpkm, imkm, jmkm, jup, jum;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double u1, u2, v1, v2, w1, w2;
	double api, aci, ami, apj, acj, amj, apk, ack, amk, mbcd, mbcu;
	double l31un, l32vn, l33wn, m31un, m32vn, m33wn, pressg, mbc=0;

	double *u,*v,*w,*p;        FLD.ptrGet(u,v,w,p);
	double *nux,*nuy,*nuz,*nu; VIS.ptrGet(nux,nuy,nuz,nu);

	# pragma omp for
	for (j=1; j<Ny; j++) { jup = (j!=Ny-1); jum = (j!=1);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		ipkm = IDX(ipa[i],j,kma[k]); jpkm = IDX(i,j+1,kma[k]);
		imkm = IDX(ima[i],j,kma[k]); jmkm = IDX(i,j-1,kma[k]);
		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

		// interpolate viscosity and velocity at step n to the position needed
		vis1 = nuy[idx]; vis2 = nuy[ip];
		vis3 = nux[idx]; vis4 = nux[jp];
		vis5 = nu[km];   vis6 = nu[idx];

		u2 = 0.5 * (u[ip] + u[ipkm]);
		u1 = 0.5 * (u[idx]+ u[km]);
		v2 = 0.5 * (v[jp] + v[jpkm]);
		v1 = 0.5 * (v[idx]+ v[km]);
		w2 = 0.5 * (w[idx]+ w[kp]);
		w1 = 0.5 * (w[idx]+ w[km]);

		// viscous terms
		api = 0.5/dx2 * vis2;
		aci =-0.5/dx2 *(vis2+vis1);
		ami = 0.5/dx2 * vis1;
		apj = 0.5/dy[j] * vis4/h[j+1];
		acj =-0.5/dy[j] *(vis4/h[j+1]+vis3/h[j]);
		amj = 0.5/dy[j] * vis3/h[j];
		apk = 1.0/dz2 * vis6;
		ack =-1.0/dz2 *(vis6+vis5);
		amk = 1.0/dz2 * vis5;

		l31un = ( vis2 * (u[ip]-u[ipkm]) - vis1 * (u[idx]-u[km]) ) / (2.0*dz*dx);
		l32vn = ( vis4 * (v[jp]-v[jpkm]) - vis3 * (v[idx]-v[km]) ) / (2.0*dz*dy[j]);
		l33wn =	api*w[ip] + aci*w[idx] + ami*w[im]
			+	apj*w[jp] + acj*w[idx] + amj*w[jm]
			+	apk*w[kp] + ack*w[idx] + amk*w[km];

		// non-linear terms
		// m33wn
		api = 0.25/dx * u2     - api;
		aci = 0.25/dx *(u2-u1) - aci;
		ami =-0.25/dx * u1     - ami;
		apj = 0.25/h[j+1]* v2                                   - apj;
		acj = 0.25/dy[j] *(v2*dy[j+1]/h[j+1] - v1*dy[j-1]/h[j]) - acj;
		amj =-0.25/h[j]  * v1                                   - amj;
		apk = 0.5/dz * w2     - apk;
		ack = 0.5/dz *(w2-w1) - ack;
		amk =-0.5/dz * w1     - amk;

		// mbcd = amj * wbc[IDX(i,0,k)];
		// mbcu = apj * wbc[IDX(i,1,k)];
		// mbc = (jum-1) * mbcd + (jup-1) * mbcu;
		apj *= jup;
		amj *= jum;
		m33wn =	api*w[ip] + aci*w[idx] + ami*w[im]
			+	apj*w[jp] + acj*w[idx] + amj*w[jm]
			+	apk*w[kp] + ack*w[idx] + amk*w[km];

		// m32vn
		w2 = ( w[idx]*dy[j+1] + w[jp]*dy[j] ) / (2.0*h[j+1]);
		w1 = ( w[idx]*dy[j-1] + w[jm]*dy[j] ) / (2.0*h[j]);

		// mbcd = .5/dy[1]    * (-w1 * (vbc[IDX(i,0,k)]+vbc[IDX(i,0,kma[k])])/2. + vis3 * (vbc[IDX(i,0,k)]-vbc[IDX(i,0,kma[k])])/dz );
		// mbcu = .5/dy[Ny-1] * ( w2 * (vbc[IDX(i,1,k)]+vbc[IDX(i,1,kma[k])])/2. - vis4 * (vbc[IDX(i,1,k)]-vbc[IDX(i,1,kma[k])])/dz );
		// mbc += (jum-1) * mbcd + (jup-1) * mbcu;
		w2 *= jup;	vis4 *= jup;
		w1 *= jum;	vis3 *= jum;
		m32vn = .5/dy[j] * ( w2*v2 - w1*v1 - ( vis4 * (v[jp]-v[jpkm]) - vis3 * (v[idx]-v[km]) )/dz ); // = (w2*v2 - w1*v1) / (2.0*dy[j]) - ( vis4 * (v[jp]-v[jpkm]) - vis3 * (v[idx]-v[km]) ) / (2.0*dz*dy[j]);

		// m31un
		w2 = 0.5 * ( w[idx] + w[ip] );
		w1 = 0.5 * ( w[idx] + w[im] );
		m31un = (w2*u2 - w1*u1) / (2.0*dx) - l31un;

		// pressure gradient term
		pressg = (p[idx] - p[km]) / dz;

		// vis3 = nux[idx];
		// vis4 = nux[jp];
		// v2 = 0.5 * (v[jp] + v[jpkm]);
		// v1 = 0.5 * (v[idx]+ v[km]);
		// amj =-0.25/h[j]  * v1 - 0.5/dy[j] * vis3/h[j];
		// apj = 0.25/h[j+1]* v2 - 0.5/dy[j] * vis4/h[j+1];
		// w2 = ( w[idx]*dy[j+1] + w[jp]*dy[j] ) / (2.0*h[j+1]);
		// w1 = ( w[idx]*dy[j-1] + w[jm]*dy[j] ) / (2.0*h[j]);
		// mbcd = amj * wbc[IDX(i,0,k)] + .5/dy[1]    * (-w1 * (vbc[IDX(i,0,k)]+vbc[IDX(i,0,kma[k])])/2. + vis3 * (vbc[IDX(i,0,k)]-vbc[IDX(i,0,kma[k])])/dz );
		// mbcu = apj * wbc[IDX(i,1,k)] + .5/dy[Ny-1] * ( w2 * (vbc[IDX(i,1,k)]+vbc[IDX(i,1,kma[k])])/2. - vis4 * (vbc[IDX(i,1,k)]-vbc[IDX(i,1,kma[k])])/dz );
		// mbc = (jum-1) * mbcd + (jup-1) * mbcu;

		// R_3 without boundary modification
		rwh[idx] = (l31un + l32vn + l33wn)
		         - (m31un + m32vn + m33wn)
		         - pressg + fbz[idx] + mbc;
	}}}
}


void IDM::mbc(Vctr &UH, const Feld &FLD, const Feld &VIS, const Feld &BC)
/* boundary modification accounting for the boundary information removed from LHS */
{
	int i, k, j0, j1, j2, jm, jn, imj0, imj1, imjn, kmj0, kmj1, kmjn;
	double u1, u2, ub, v1, v2, vb, w1, w2, wb, vis3, vis4, amj, apj;

	double *ruh,*rvh,*rwh;      UH.ptrGet(ruh,rvh,rwh);
	double *u,*v,*w,*p;        FLD.ptrGet(u,v,w,p);
	double *nux,*nuy,*nuz,*nu; VIS.ptrGet(nux,nuy,nuz,nu);
	double *ubc,*vbc,*wbc,*pbc; BC.ptrGet(ubc,vbc,wbc,pbc);

	# pragma omp for
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		j0 = IDX(i,0 ,k); imj0 = IDX(ima[i],0 ,k); kmj0 = IDX(i,0 ,kma[k]);
		j1 = IDX(i,1 ,k); imj1 = IDX(ima[i],1 ,k); kmj1 = IDX(i,1 ,kma[k]);
		j2 = IDX(i,2 ,k);
		jm = IDX(i,Ny-1,k);
		jn = IDX(i,Ny,k); imjn = IDX(ima[i],Ny,k); kmjn = IDX(i,Ny,kma[k]);

		// mbc_1, j = 1
		vis3 = nuz[j1];
		v1 = .5 * (v[j1] + v[imj1]);	amj = -.25/h[1] * v1 - .5/dy[1] * vis3/h[1];
		u1 = .5/h[1] * (u[j1]*dy[0] + u[j0]*dy[1]);
		ub = ubc[j0];
		vb = .5 * (vbc[j0] + vbc[imj0]);
		ruh[j1] -= amj * ub + .5/dy[1] * ( -u1 * vb + vis3 * (vbc[j0]-vbc[imj0])/dx );
		// mbc_1, j = Ny-1
		vis4 = nuz[jn];
		v2 = .5 * (v[jn] + v[imjn]);	apj = .25/h[Ny]* v2 - .5/dy[Ny-1] * vis4/h[Ny];
		u2 = .5/h[Ny] * (u[jm]*dy[Ny] + u[jn]*dy[Ny-1]);
		ub = ubc[j1];
		vb = .5 * (vbc[j1] + vbc[imj1]);
		ruh[jm] -= apj * ub + .5/dy[Ny-1] * ( u2 * vb - vis4 * (vbc[j1]-vbc[imj1])/dx );

		// mbc_2, j = 2
		vis3 = nu[j1];
		v1 = .5 * (v[j2] + v[j1]);	amj = -.5/h[2] * v1 - 1./h[2] * vis3/dy[1];
		rvh[j2] -= amj * vbc[j0];
		// mbc_2, j = Ny-1
		vis4 = nu[jm];
		v2 = .5 * (v[jm] + v[jn]);	apj = .5/h[Ny-1] * v2 - 1./h[Ny-1] * vis4/dy[Ny-1];
		rvh[jm] -= apj * vbc[j1];

		// mbc_3, j = 1
		vis3 = nux[j1];
		v1 = .5 * (v[j1] + v[kmj1]);	amj = -.25/h[1] * v1 - .5/dy[1] * vis3/h[1];
		w1 = .5/h[1] * (w[j1]*dy[0] + w[j0]*dy[1]);
		wb = wbc[j0];
		vb = .5 * (vbc[j0]+vbc[kmj0]);
		rwh[j1] -= amj * wb + .5/dy[1] * ( - w1*vb + vis3 * (vbc[j0]-vbc[kmj0])/dz );
		// mbc_3, j = Ny-1
		vis4 = nux[jn];
		v2 = .5 * (v[jn] + v[kmjn]);	apj = .25/h[Ny] * v2 - .5/dy[Ny-1] * vis4/h[Ny];
		w2 = .5/h[Ny] * (w[jm]*dy[Ny] + w[jn]*dy[Ny-1]);
		wb = wbc[j1];
		vb = .5 * (vbc[j1] + vbc[kmj1]);
		rwh[jm] -= apj * wb + .5/dy[Ny-1] * ( w2*vb - vis4 * (vbc[j1]-vbc[kmj1])/dz );
	}}
}






void IDM::getuh1(Vctr &UH, const Vctr &U, const Feld &VIS, double dt)
/* compute deltaU^**, result returned by uh (the RHS ruh should be pre stored in uh ) */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, imjp, imkp, imjm, imkm;
	double u1, u2, v1, v2, w1, w2;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double *api = new double [Nx], *aci = new double [Nx], *ami = new double [Nx], *R1 = new double [Nx];
	double *apj = new double [Ny], *acj = new double [Ny], *amj = new double [Ny], *R2 = new double [Ny];
	double *apk = new double [Nz], *ack = new double [Nz], *amk = new double [Nz], *R3 = new double [Nz];
	Matrix matx(Nx);
	Matrix maty(Ny-1);
	Matrix matz(Nz);

	double *uh,*vh,*wh;         UH.ptrGet(uh,vh,wh);
	double *u,*v,*w;             U.ptrGet(u,v,w);
	double *nux,*nuy,*nuz,*nu; VIS.ptrGet(nux,nuy,nuz,nu);

	# pragma omp for
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		for (j=1; j<Ny; j++) {
			idx = IDX(i,j,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
			imjp = IDX(ima[i],j+1,k); imjm = IDX(ima[i],j-1,k);

			v2 = 0.5 * ( v[jp] + v[imjp] );
			v1 = 0.5 * ( v[idx] + v[im] );
			vis3 = nuz[idx];
			vis4 = nuz[jp];

			apj[j] = ( 0.25/h[j+1]* v2                                   - 0.5/dy[j] * vis4/h[j+1]            ) * dt;
			acj[j] = ( 0.25/dy[j] *(v2*dy[j+1]/h[j+1] - v1*dy[j-1]/h[j]) + 0.5/dy[j] *(vis4/h[j+1]+vis3/h[j]) ) * dt + 1;
			amj[j] = (-0.25/h[j]  * v1                                   - 0.5/dy[j] * vis3/h[j]              ) * dt;

			R2 [j] = dt * uh[idx];
		}
		maty.tdma( & amj[1], & acj[1], & apj[1], & R2[1] ); // apj at j=Ny-1 and amj at j=1 are redundant in tdma, thus no need for explicit removal
		for (j=1; j<Ny; j++) uh[IDX(i,j,k)] = R2[j];
	}}

	# pragma omp for
	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			idx = IDX(i,j,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

			u2 = 0.5 * ( u[idx] + u[ip] );
			u1 = 0.5 * ( u[idx] + u[im] );
			vis1 = nu[im];
			vis2 = nu[idx];

			api[i] = ( 0.5/dx * u2     - 1.0/dx2 * vis2       ) * dt;
			aci[i] = ( 0.5/dx *(u2-u1) + 1.0/dx2 *(vis2+vis1) ) * dt + 1;
			ami[i] = (-0.5/dx * u1     - 1.0/dx2 * vis1       ) * dt;

			R1 [i] = uh[idx];
		}
		matx.ctdma( ami, aci, api, R1 );
		for (i=0; i<Nx; i++) uh[IDX(i,j,k)] = R1[i];
	}}

	# pragma omp for
	for (j=1; j<Ny; j++) {
	for (i=0; i<Nx; i++) {
		for (k=0; k<Nz; k++) {
			idx = IDX(i,j,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
			imkp = IDX(ima[i],j,kpa[k]); imkm = IDX(ima[i],j,kma[k]);

			w2 = 0.5 * ( w[kp] + w[imkp] );
			w1 = 0.5 * ( w[idx] + w[im] );
			vis5 = nuy[idx];
			vis6 = nuy[kp];

			apk[k] = ( 0.25/dz * w2     - 0.5/dz2 * vis6       ) * dt;
			ack[k] = ( 0.25/dz *(w2-w1) + 0.5/dz2 *(vis6+vis5) ) * dt + 1;
			amk[k] = (-0.25/dz * w1     - 0.5/dz2 * vis5       ) * dt;

			R3 [k] = uh[idx];
		}
		matz.ctdma( amk, ack, apk, R3 );
		for (k=0; k<Nz; k++) uh[IDX(i,j,k)] = R3[k];
	}}

	delete [] api; delete [] aci; delete [] ami; delete [] R1;
	delete [] apj; delete [] acj; delete [] amj; delete [] R2;
	delete [] apk; delete [] ack; delete [] amk; delete [] R3;
}


void IDM::getuh2(Vctr &UH, const Vctr &U, const Feld &VIS, double dt)
/* compute deltaU^**, result returned by vh (the RHS rvh should be pre stored in uh ) */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, ipjm, jmkp, imjm, jmkm;
	double u1, u2, v1, v2, w1, w2, l21uh, m21uh;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double *api = new double [Nx], *aci = new double [Nx], *ami = new double [Nx], *R1 = new double [Nx];
	double *apj = new double [Ny], *acj = new double [Ny], *amj = new double [Ny], *R2 = new double [Ny];
	double *apk = new double [Nz], *ack = new double [Nz], *amk = new double [Nz], *R3 = new double [Nz];
	Matrix matx(Nx);
	Matrix maty(Ny-2);
	Matrix matz(Nz);

	double *uh,*vh,*wh;         UH.ptrGet(uh,vh,wh);
	double *u,*v,*w;             U.ptrGet(u,v,w);
	double *nux,*nuy,*nuz,*nu; VIS.ptrGet(nux,nuy,nuz,nu);

	# pragma omp for
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		for (j=2; j<Ny; j++) {
			idx = IDX(i,j,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
			ipjm = IDX(ipa[i],j-1,k); imjm = IDX(ima[i],j-1,k);
		
			v2 = 0.5 * (v[idx] + v[jp]);
			v1 = 0.5 * (v[idx] + v[jm]);
			vis3 = nu[jm];
			vis4 = nu[idx];

			apj[j] = ( 0.5/h[j]* v2     - 1.0/h[j] * vis4/dy[j]               ) * dt;
			acj[j] = ( 0.5/h[j]*(v2-v1) + 1.0/h[j] *(vis4/dy[j]+vis3/dy[j-1]) ) * dt + 1;
			amj[j] = (-0.5/h[j]* v1     - 1.0/h[j] * vis3/dy[j-1]             ) * dt;

			// m21uh
			u2 = ( uh[ip]*dy[j-1] + uh[ipjm]*dy[j] ) / (2.0*h[j]);
			u1 = ( uh[idx]*dy[j-1] + uh[jm]*dy[j] ) / (2.0*h[j]);
			v2 = 0.5 * ( v[idx] + v[ip] );
			v1 = 0.5 * ( v[idx] + v[im] );
			vis1 = 0.25 * ( (nu[idx]+nu[im]) * dy[j-1] + (nu[jm]+nu[imjm]) * dy[j] ) / h[j];
			vis2 = 0.25 * ( (nu[idx]+nu[ip]) * dy[j-1] + (nu[jm]+nu[ipjm]) * dy[j] ) / h[j];
			
			l21uh = ( vis2 * (uh[ip]-uh[ipjm]) - vis1 * (uh[idx]-uh[jm]) ) / (2.0*h[j]*dx);
			m21uh = (v2*u2 - v1*u1) / (2.0*dx) - l21uh;

			R2 [j] = dt * ( vh[idx] - m21uh );
		}
		maty.tdma( & amj[2], & acj[2], & apj[2], & R2[2] ); // apj at j=Ny-1 and amj at j=1 are redundant in tdma, thus no need for explicit removal
		for (j=2; j<Ny; j++) vh[IDX(i,j,k)] = R2[j];
	}}

	# pragma omp for
	for (j=2; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			idx = IDX(i,j,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
			ipjm = IDX(ipa[i],j-1,k); imjm = IDX(ima[i],j-1,k);

			u2 = ( u[ip]*dy[j-1] + u[ipjm]*dy[j] ) / (2.0*h[j]);
			u1 = ( u[idx]*dy[j-1] + u[jm]*dy[j] ) / (2.0*h[j]);
			vis1 = nuz[idx];
			vis2 = nuz[ip];

			api[i] = ( 0.25/dx * u2     - 0.5/dx2 * vis2       ) * dt;
			aci[i] = ( 0.25/dx *(u2-u1) + 0.5/dx2 *(vis2+vis1) ) * dt + 1;
			ami[i] = (-0.25/dx * u1     - 0.5/dx2 * vis1       ) * dt;

			R1 [i] = vh[idx];
		}
		matx.ctdma( ami, aci, api, R1 );
		for (i=0; i<Nx; i++) vh[IDX(i,j,k)] = R1[i];
	}}

	# pragma omp for
	for (j=2; j<Ny; j++) {
	for (i=0; i<Nx; i++) {
		for (k=0; k<Nz; k++) {
			idx = IDX(i,j,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
			jmkp = IDX(i,j-1,kpa[k]); jmkm = IDX(i,j-1,kma[k]);

			w2 = ( w[kp]*dy[j-1] + w[jmkp]*dy[j] ) / (2.0*h[j]);
			w1 = ( w[idx]*dy[j-1] + w[jm]*dy[j] ) / (2.0*h[j]);
			vis5 = nux[idx];
			vis6 = nux[kp];

			apk[k] = ( 0.25/dz * w2     - 0.5/dz2 * vis6       ) * dt;
			ack[k] = ( 0.25/dz *(w2-w1) + 0.5/dz2 *(vis6+vis5) ) * dt + 1;
			amk[k] = (-0.25/dz * w1     - 0.5/dz2 * vis5       ) * dt;

			R3 [k] = vh[idx];
		}
		matz.ctdma( amk, ack, apk, R3 );
		for (k=0; k<Nz; k++) vh[IDX(i,j,k)] = R3[k];
	}}
	
	delete [] api; delete [] aci; delete [] ami; delete [] R1;
	delete [] apj; delete [] acj; delete [] amj; delete [] R2;
	delete [] apk; delete [] ack; delete [] amk; delete [] R3;
}


void IDM::getuh3(Vctr &UH, const Vctr &U, const Feld &VIS, double dt)
/* compute deltaW^*, deltaV^*, deltaU^*, and update U^*, V^*, W^* (the RHS rwh should be pre stored in wh ) */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, jup, jum;
	int ipkm, jpkm, imjp, imkp, jmkp, imkm, jmkm, imjm;
	double u1, u2, v1, v2, w1, w2;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double l31uh, l32vh, l23wh, l12vh, l13wh;
	double m31uh, m32vh, m23wh, m12vh, m13wh;
	double *api = new double [Nx], *aci = new double [Nx], *ami = new double [Nx], *R1 = new double [Nx];
	double *apj = new double [Ny], *acj = new double [Ny], *amj = new double [Ny], *R2 = new double [Ny];
	double *apk = new double [Nz], *ack = new double [Nz], *amk = new double [Nz], *R3 = new double [Nz];
	Matrix matx(Nx);
	Matrix maty(Ny-1);
	Matrix matz(Nz);

	double *uh,*vh,*wh;         UH.ptrGet(uh,vh,wh);
	double *u,*v,*w;             U.ptrGet(u,v,w);
	double *nux,*nuy,*nuz,*nu; VIS.ptrGet(nux,nuy,nuz,nu);

	// ( I + dt M_33^2 )
	# pragma omp for
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		for (j=1; j<Ny; j++) {  jup = (j!=Ny-1); jum = (j!=1);
			idx = IDX(i,j,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
			ipkm = IDX(ipa[i],j,kma[k]); jpkm = IDX(i,j+1,kma[k]);
			imkm = IDX(ima[i],j,kma[k]); jmkm = IDX(i,j-1,kma[k]);

			v2 = 0.5 * ( v[jp] + v[jpkm] );
			v1 = 0.5 * ( v[idx] + v[km] );
			vis3 = nux[idx];
			vis4 = nux[jp];

			apj[j] = ( 0.25/h[j+1]* v2                                   - 0.5/dy[j] * vis4/h[j+1]            ) * dt;
			acj[j] = ( 0.25/dy[j] *(v2*dy[j+1]/h[j+1] - v1*dy[j-1]/h[j]) + 0.5/dy[j] *(vis4/h[j+1]+vis3/h[j]) ) * dt + 1;
			amj[j] = (-0.25/h[j]  * v1                                   - 0.5/dy[j] * vis3/h[j]              ) * dt;

			// m31uh
			u2 = 0.5 * ( uh[ip] + uh[ipkm] );
			u1 = 0.5 * ( uh[idx] + uh[km] );
			w2 = 0.5 * ( w[idx] + w[ip] );
			w1 = 0.5 * ( w[idx] + w[im] );
			vis1 = 0.25 * ( nu[idx] + nu[im] + nu[km] + nu[imkm] );
			vis2 = 0.25 * ( nu[idx] + nu[ip] + nu[km] + nu[ipkm] );

			l31uh = ( vis2 * (uh[ip]-uh[ipkm]) - vis1 * (uh[idx]-uh[km]) ) / (2.0*dz*dx);
			m31uh = (w2*u2 - w1*u1) / (2.0*dx) - l31uh;

			// m32vh
			v2 = 0.5 * ( vh[jp] + vh[jpkm] );
			v1 = 0.5 * ( vh[idx] + vh[km] );
			w2 = ( w[idx]*dy[j+1] + w[jp]*dy[j] ) / (2.0*h[j+1]);
			w1 = ( w[idx]*dy[j-1] + w[jm]*dy[j] ) / (2.0*h[j]);

			w2 *= jup;	vis4 *= jup;
			w1 *= jum;	vis3 *= jum;
			l32vh = ( vis4 * (vh[jp]-vh[jpkm]) - vis3 * (vh[idx]-vh[km]) ) / (2.0*dz*dy[j]);
			m32vh = (w2*v2 - w1*v1) / (2.0*dy[j]) - l32vh;

			R2 [j] = dt * ( wh[idx] - m31uh - m32vh );
		}
		maty.tdma( & amj[1], & acj[1], & apj[1], & R2[1] );
		for (j=1; j<Ny; j++) wh[IDX(i,j,k)] = R2[j];
	}}

	// ( I + dt M_33^1 )
	# pragma omp for
	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			idx = IDX(i,j,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
			ipkm = IDX(ipa[i],j,kma[k]); imkm = IDX(ima[i],j,kma[k]);

			u2 = 0.5 * ( u[ip] + u[ipkm] );
			u1 = 0.5 * ( u[idx] + u[km] );
			vis1 = nuy[idx];
			vis2 = nuy[ip];

			api[i] = ( 0.25/dx * u2     - 0.5/dx2 * vis2       ) * dt;
			aci[i] = ( 0.25/dx *(u2-u1) + 0.5/dx2 *(vis2+vis1) ) * dt + 1;
			ami[i] = (-0.25/dx * u1     - 0.5/dx2 * vis1       ) * dt;

			R1 [i] = wh[idx];
		}
		matx.ctdma( ami, aci, api, R1 );
		for (i=0; i<Nx; i++) wh[IDX(i,j,k)] = R1[i];
	}}

	// ( I + dt M_33^3 )
	# pragma omp for
	for (j=1; j<Ny; j++) {
	for (i=0; i<Nx; i++) {
		for (k=0; k<Nz; k++) {
			idx = IDX(i,j,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

			w2 = 0.5 * ( w[idx] + w[kp] );
			w1 = 0.5 * ( w[idx] + w[km] );
			vis5 = nu[km];
			vis6 = nu[idx];

			apk[k] = ( 0.5/dz * w2     - 1.0/dz2 * vis6       ) * dt;
			ack[k] = ( 0.5/dz *(w2-w1) + 1.0/dz2 *(vis6+vis5) ) * dt + 1;
			amk[k] = (-0.5/dz * w1     - 1.0/dz2 * vis5       ) * dt;

			R3 [k] = wh[idx];
		}
		matz.ctdma( amk, ack, apk, R3 );
		for (k=0; k<Nz; k++) wh[IDX(i,j,k)] = R3[k];
	}}

	delete [] api; delete [] aci; delete [] ami; delete [] R1;
	delete [] apj; delete [] acj; delete [] amj; delete [] R2;
	delete [] apk; delete [] ack; delete [] amk; delete [] R3;

	// update dvh
	# pragma omp for
	for (j=2; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
		jmkp = IDX(i,j-1,kpa[k]); jmkm = IDX(i,j-1,kma[k]);

		// m23wh
		v2 = 0.5 * ( v[idx] + v[kp] );
		v1 = 0.5 * ( v[idx] + v[km] );
		w2 = ( wh[kp]*dy[j-1] + wh[jmkp]*dy[j] ) / (2.0*h[j]);
		w1 = ( wh[idx]*dy[j-1] + wh[jm]*dy[j] ) / (2.0*h[j]);
		vis5 = nux[idx];
		vis6 = nux[kp];

		l23wh = ( vis6 * (wh[kp]-wh[jmkp]) - vis5 * (wh[idx]-wh[jm]) ) / (2.0*h[j]*dz);
		m23wh = (v2*w2 - v1*w1) / (2.0*dz) - l23wh;

		vh[idx] -= dt * m23wh;
	}}}
	
	// update duh
	# pragma omp for
	for (j=1; j<Ny; j++) {  jup = (j!=Ny-1); jum = (j!=1);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
		imjp = IDX(ima[i],j+1,k); imkp = IDX(ima[i],j,kpa[k]);
		imjm = IDX(ima[i],j-1,k); imkm = IDX(ima[i],j,kma[k]);

		// m12vh
		u2 = ( u[idx]*dy[j+1] + u[jp]*dy[j] ) / (2.0*h[j+1]);
		u1 = ( u[idx]*dy[j-1] + u[jm]*dy[j] ) / (2.0*h[j]);
		v2 = 0.5 * ( vh[jp] + vh[imjp] );
		v1 = 0.5 * ( vh[idx] + vh[im] );
		vis3 = nuz[idx];
		vis4 = nuz[jp];

		u2 *= jup;	vis4 *= jup;
		u1 *= jum;	vis3 *= jum;
		l12vh = ( vis4 * (vh[jp]-vh[imjp]) - vis3 * (vh[idx]-vh[im]) ) / (2.0*dx*dy[j]);
		m12vh = (u2*v2 - u1*v1) / (2.0*dy[j]) - l12vh;

		// m13wh
		u2 = 0.5 * ( u[idx] + u[kp] );
		u1 = 0.5 * ( u[idx] + u[km] );
		w2 = 0.5 * ( wh[kp] + wh[imkp] );
		w1 = 0.5 * ( wh[idx] + wh[im] );
		vis5 = nuy[idx];
		vis6 = nuy[kp];

		l13wh = ( vis6 * (wh[kp]-wh[imkp]) - vis5 * (wh[idx]-wh[im]) ) / (2.0*dx*dz);
		m13wh = (u2*w2 - u1*w1) / (2.0*dz) - l13wh;

		uh[idx] -= dt * ( m12vh + m13wh );
	}}}

	// update intermediate velocity field
	# pragma omp for
	for (j=1; j<Ny; j++) {	jum = (j!=1);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		uh[idx] += u[idx];
		vh[idx] += v[idx] * jum; // j=1 is redundant
		wh[idx] += w[idx];
	}}}
}



/***** projector computation *****/

void IDM::rhsdp(double *rdp, const Vctr &UH, const Vctr &UBC, double dt)
/* compute RHS of Poisson equation and store in dp */
{
	int i, j, k, idx, ip, jp, kp, j1, j0, jup, jum;

	double *uh,*vh,*wh;     UH.ptrGet(uh,vh,wh);
	double *ubc,*vbc,*wbc; UBC.ptrGet(ubc,vbc,wbc);

	for (j=1; j<Ny; j++) { jup = (j!=Ny-1); jum = (j!=1);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx= IDX(i,j,k); ip = IDX(ipa[i],j,k);
		j1 = IDX(i,1,k); jp = IDX(i,j+1,k);
		j0 = IDX(i,0,k); kp = IDX(i,j,kpa[k]);

		// ( Du^h - cbc ) / dt
		rdp[idx] = 1.0/dt * (
			( uh[ip] - uh[idx] ) / dx
		+	( wh[kp] - wh[idx] ) / dz
		+	( vh[jp] * jup + vbc[j1] * (1-jup)
			- vh[idx]* jum - vbc[j0] * (1-jum) ) / dy[j]	);
	}}}
}

void IDM::getfdp(double *fdp, double refp)
/* compute FDP (in Fourier space), result returned by fdp (the RHS frdp should be pre stored in fdp ) */
{
	int i, j, k, idx;
	int Nxc = (int) (Nx/2+1), Nxr = 2 * Nxc, Nxzr = Nz * Nxr;
	double *cpj = new double [Ny];
	double *ccj = new double [Ny];
	double *cmj = new double [Ny];
	double *cfj1= new double [Ny], *cfj2= new double [Ny];
	Matrix matyp(Ny-1);

	for (k=0; k<Nz; k++) {
	for (i=0; i<Nxc; i++) { // negative k_x need not solve
		for (j=1; j<Ny; j++) {
			idx = Nxzr * j + Nxr * k + 2*i;

			cpj[j] = ppj[j];
			ccj[j] = pcj[j] - ak3[k] - ak1[i];
			cmj[j] = pmj[j];
			cfj1[j] = fdp[idx];		// real part
			cfj2[j] = fdp[idx+1];	// imaginary part
		}

		/* note:
			In Huang's program, imag part of FDP at wavenumbers
			(kx=0,kz=0), (kx=0,kz=Nz/2), (kx=Nx/2,kz=0), (kx=Nx/2,kz=Nz/2)
			are explicitly set to 0. This is not needed here, since
			zero imag FRDP lead directly to zero imag FDP	*/

		// set reference pressure P(kx=0,kz=0,j=1) to 0
		if (k==0 && i==0) {
			cpj[1] = 0;
			ccj[1] = 1;
			cmj[1] = 0;
			cfj1[1] = - Nxz * refp;	// fdp(kx=0,kz=0,j=1) = -Nxz * fp(kx=0,kz=0,j=1) ==> <dp(j=1)> = - <p(j=1)>
			cfj2[1] = 0;
		}

		matyp.tdma( & cmj[1], & ccj[1], & cpj[1], & cfj1[1] );
		matyp.tdma( & cmj[1], & ccj[1], & cpj[1], & cfj2[1] );

		for (j=1; j<Ny; j++) {
			idx = Nxzr * j + Nxr * k + 2*i;
			fdp[idx] = cfj1[j];
			fdp[idx+1] = cfj2[j];
		}
	}}

	delete [] cpj; delete [] ccj; delete [] cmj;
	delete [] cfj1;delete [] cfj2;
}











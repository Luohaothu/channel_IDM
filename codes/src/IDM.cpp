# include <iostream>
# include <math.h>
# include <omp.h>

# include "IDM.h"
# include "Matrix.h"

using namespace std;



/***** configuration *****/

int IDM::ompset(int n)
{
	int nprocs = omp_get_num_procs();
	if (n > 0) this->nthrds = (n > nprocs ? nprocs : n);
	else {
	/* automatically decide number of OpenMP threads based on grid number */
		int ngrids = log2(Nx*Ny*Nz) - 16;
		nprocs = (nprocs > 32 ? 32 : nprocs);
		ngrids = (ngrids < 1 ? 1 : ngrids);
		this->nthrds = (ngrids > nprocs ? nprocs : ngrids);
		// cout << "\nnumber of OpenMP threads for uhcalc: " << nthrds << endl;
	}
	return this->nthrds;
}


/***** computation interfaces *****/

void IDM::uhcalc(Vctr &UH, Vctr &U, Scla &P, Vctr &UBC, Scla &NU, double mpg[3])
{
	double *u   =   U.bulkGet(1), *v   =   U.bulkGet(2), *w   =   U.bulkGet(3), *p  =  P.bulkGet();
	double *uh  =  UH.bulkGet(1), *vh  =  UH.bulkGet(2), *wh  =  UH.bulkGet(3);
	double *ubc = UBC.bulkGet(1), *vbc = UBC.bulkGet(2), *wbc = UBC.bulkGet(3), *nu = NU.bulkGet();

	omp_set_num_threads(this->nthrds);

	# pragma omp parallel
	{
		// calculate RUH (which shares memory with UH)
		urhs1(uh, u, v, w, p, nu, mpg[0]);
		urhs2(vh, u, v, w, p, nu        );
		urhs3(wh, u, v, w, p, nu, mpg[2]);
		mbc(uh, vh, wh, u, v, w, ubc, vbc, wbc, nu);
		// calculate UH
		getuh1(uh,         u, v, w, nu);
		getuh2(uh, vh,     u, v, w, nu);
		getuh3(uh, vh, wh, u, v, w, nu);
	}
}

void IDM::dpcalc(Scla &DP, Vctr &UH, Scla &P, Vctr &UBC)
{
	double *uh = UH.bulkGet(1), *vh = UH.bulkGet(2), *wh = UH.bulkGet(3);
	double *dp = DP.bulkGet(), *fdp = DP.bulkGetF(), *vbc = UBC.bulkGet(2);
	double refp = P.layerMean(1);

	rhsdp(dp, uh, vh, wh, vbc);	// rdp (which shares memory with dp)
	DP.fft();					// rdp->frdp
	getfdp(fdp, refp);			// frdp->fdp
	DP.ifft();					// fdp->dp
}

void IDM::upcalc(Vctr &U, Scla &P, Vctr &UH, Scla &DP, double *mpg)
{
	double *u  =  U.bulkGet(1), *v  =  U.bulkGet(2), *w  =  U.bulkGet(3), *p  =  P.bulkGet();
	double *uh = UH.bulkGet(1), *vh = UH.bulkGet(2), *wh = UH.bulkGet(3), *dp = DP.bulkGet();

	update(u, v, w, p, uh, vh, wh, dp);
	meanpg(mpg, U);
}



/***** computations functions *****/

void IDM::urhs1(double *ruh, double *u, double *v, double *w, double *p, double *nu, double mpg1)
/* compute right hand side R_1 for intermediate velocity at all non-wall grid points */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, imjp, imkp, imjm, imkm, jup, jum;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double u1, u2, v1, v2, w1, w2;
	double api, aci, ami, apj, acj, amj, apk, ack, amk;
	double l11un, l12vn, l13wn, m11un, m12vn, m13wn, pressg;

	# pragma omp for
	for (j=1; j<Ny; j++) {	jup = j==Ny-1 ? 0 : 1;	jum = j==1 ? 0 : 1; // indicate the secondary boundary
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		imjp = IDX(ima[i],j+1,k); imkp = IDX(ima[i],j,kpa[k]);
		imjm = IDX(ima[i],j-1,k); imkm = IDX(ima[i],j,kma[k]);
		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

		// interpolate viscosity and velocity at step n to the position needed
		vis1 = nu[im];
		vis2 = nu[idx];
		vis3 = 0.25 * ( (nu[idx]+nu[im]) * dy[j-1] + (nu[jm]+nu[imjm]) * dy[j] ) / h[j];
		vis4 = 0.25 * ( (nu[idx]+nu[im]) * dy[j+1] + (nu[jp]+nu[imjp]) * dy[j] ) / h[j+1];
		vis5 = 0.25 * ( nu[idx] + nu[im] + nu[km] + nu[imkm] );
		vis6 = 0.25 * ( nu[idx] + nu[im] + nu[kp] + nu[imkp] );

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

		apj *= jup;
		amj *= jum;
		m11un =	api*u[ip] + aci*u[idx] + ami*u[im]
			+	apj*u[jp] + acj*u[idx] + amj*u[jm]
			+	apk*u[kp] + ack*u[idx] + amk*u[km];

		// m12vn
		u2 = ( u[idx]*dy[j+1] + u[jp]*dy[j] ) / (2.0*h[j+1]);
		u1 = ( u[idx]*dy[j-1] + u[jm]*dy[j] ) / (2.0*h[j]);

		u2 *= jup;	vis4 *= jup;
		u1 *= jum;	vis3 *= jum;
		m12vn = (u2*v2 - u1*v1) / (2.0*dy[j]) - ( vis4 * (v[jp]-v[imjp]) - vis3 * (v[idx]-v[im]) ) / (2.0*dx*dy[j]);

		// m13wn
		u2 = 0.5 * ( u[idx] + u[kp] );
		u1 = 0.5 * ( u[idx] + u[km] );
		m13wn = (u2*w2 - u1*w1) / (2.0*dz) - l13wn;

		// pressure gradient term
		pressg = (p[idx] - p[im]) / dx + mpg1;

		// R_1 without boundary modification
		ruh[idx] = (l11un + l12vn + l13wn) - (m11un + m12vn + m13wn) - pressg;
	}}}
}

void IDM::urhs2(double *rvh, double *u, double *v, double *w, double *p, double *nu)
/* compute right hand side R_2 for intermediate velocity at all non-wall grid points */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, ipjm, jmkp, imjm, jmkm, jup, jum;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double u1, u2, v1, v2, w1, w2;
	double api, aci, ami, apj, acj, amj, apk, ack, amk;
	double l21un, l22vn, l23wn, m21un, m22vn, m23wn, pressg;

	# pragma omp for
	for (j=2; j<Ny; j++) {	jup = j==Ny-1 ? 0 : 1;	jum = j==2 ? 0 : 1;
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		ipjm = IDX(ipa[i],j-1,k); jmkp = IDX(i,j-1,kpa[k]);
		imjm = IDX(ima[i],j-1,k); jmkm = IDX(i,j-1,kma[k]);
		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

		// interpolate viscosity and velocity at step n to the position needed
		vis1 = 0.25 * ( (nu[idx]+nu[im]) * dy[j-1] + (nu[jm]+nu[imjm]) * dy[j] ) / h[j];
		vis2 = 0.25 * ( (nu[idx]+nu[ip]) * dy[j-1] + (nu[jm]+nu[ipjm]) * dy[j] ) / h[j];
		vis3 = nu[jm];
		vis4 = nu[idx];
		vis5 = 0.25 * ( (nu[idx]+nu[km]) * dy[j-1] + (nu[jm]+nu[jmkm]) * dy[j] ) / h[j];
		vis6 = 0.25 * ( (nu[idx]+nu[kp]) * dy[j-1] + (nu[jm]+nu[jmkp]) * dy[j] ) / h[j];

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

		// R_2 without boundary modification
		rvh[idx] = (l21un + l22vn + l23wn) - (m21un + m22vn + m23wn) - pressg;
	}}}
}

void IDM::urhs3(double *rwh, double *u, double *v, double *w, double *p, double *nu, double mpg3)
/* compute right hand side R_3 for intermediate velocity at all non-wall grid points */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, ipkm, jpkm, imkm, jmkm, jup, jum;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double u1, u2, v1, v2, w1, w2;
	double api, aci, ami, apj, acj, amj, apk, ack, amk;
	double l31un, l32vn, l33wn, m31un, m32vn, m33wn, pressg;

	# pragma omp for
	for (j=1; j<Ny; j++) {	jup = j==Ny-1 ? 0 : 1;	jum = j==1 ? 0 : 1;
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		ipkm = IDX(ipa[i],j,kma[k]); jpkm = IDX(i,j+1,kma[k]);
		imkm = IDX(ima[i],j,kma[k]); jmkm = IDX(i,j-1,kma[k]);
		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

		// interpolate viscosity and velocity at step n to the position needed
		vis1 = 0.25 * ( nu[idx] + nu[im] + nu[km] + nu[imkm] );
		vis2 = 0.25 * ( nu[idx] + nu[ip] + nu[km] + nu[ipkm] );
		vis3 = 0.25 * ( (nu[idx]+nu[km]) * dy[j-1] + (nu[jm]+nu[jmkm]) * dy[j] ) / h[j];
		vis4 = 0.25 * ( (nu[idx]+nu[km]) * dy[j+1] + (nu[jp]+nu[jpkm]) * dy[j] ) / h[j+1];
		vis5 = nu[km];
		vis6 = nu[idx];

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

		apj *= jup;
		amj *= jum;
		m33wn =	api*w[ip] + aci*w[idx] + ami*w[im]
			+	apj*w[jp] + acj*w[idx] + amj*w[jm]
			+	apk*w[kp] + ack*w[idx] + amk*w[km];

		// m32vn
		w2 = ( w[idx]*dy[j+1] + w[jp]*dy[j] ) / (2.0*h[j+1]);
		w1 = ( w[idx]*dy[j-1] + w[jm]*dy[j] ) / (2.0*h[j]);

		w2 *= jup;	vis4 *= jup;
		w1 *= jum;	vis3 *= jum;
		m32vn = (w2*v2 - w1*v1) / (2.0*dy[j]) - ( vis4 * (v[jp]-v[jpkm]) - vis3 * (v[idx]-v[km]) ) / (2.0*dz*dy[j]);

		// m31un
		w2 = 0.5 * ( w[idx] + w[ip] );
		w1 = 0.5 * ( w[idx] + w[im] );
		m31un = (w2*u2 - w1*u1) / (2.0*dx) - l31un;

		// pressure gradient term
		pressg = (p[idx] - p[km]) / dz + mpg3;

		// R_3 without boundary modification
		rwh[idx] = (l31un + l32vn + l33wn) - (m31un + m32vn + m33wn) - pressg;
	}}}
}


void IDM::mbc(
	double *ruh,double *rvh,double *rwh,
	double *u,	double *v,	double *w,
	double *ubc,double *vbc,double *wbc, double *nu	)
/* boundary modification accounting for the boundary information removed from LHS */
{
	int i, k, j0, j1, j2, jm, jn, imj0, imj1, imjn, kmj0, kmj1, kmjn;
	double u1, u2, ub, v1, v2, vb, w1, w2, wb, vis3, vis4;

	# pragma omp for
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		j0 = IDX(i,0 ,k); imj0 = IDX(ima[i],0 ,k); kmj0 = IDX(i,0 ,kma[k]);
		j1 = IDX(i,1 ,k); imj1 = IDX(ima[i],1 ,k); kmj1 = IDX(i,1 ,kma[k]);
		j2 = IDX(i,2 ,k);
		jm = IDX(i,Ny-1,k);
		jn = IDX(i,Ny,k); imjn = IDX(ima[i],Ny,k); kmjn = IDX(i,Ny,kma[k]);
	
		// mbc_1, j = Ny-1
		u2 = u[jn];   v2 = 0.5 * (v[jn] + v[imjn]);
		ub = ubc[j1]; vb = 0.5 * (vbc[j1] + vbc[imj1]);
		vis4 = 0.5 * (nu[jn] + nu[imjn]);
		ruh[jm] -= 0.5/dy[Ny-1] * ( v2*ub + u2*vb - vis4*ub/h[Ny] - vis4*(v[jn]-v[imjn])/dx );
		// mbc_1, j = 1
		u1 = u[j0];   v1 = 0.5 * (v[j1] + v[imj1]);
		ub = ubc[j0]; vb = 0.5 * (vbc[j0] + vbc[imj0]);
		vis3 = 0.5 * (nu[j0] + nu[imj0]);
		ruh[j1] -= 0.5/dy[1] * ( - v1*ub - u1*vb - vis3*ub/h[1] + vis3*(v[j1]-v[imj1])/dx );

		// mbc_2, j = Ny-1
		v2 = 0.5 * (v[jn] + v[jm]);
		vb = vbc[j1];
		vis4 = nu[jn];
		rvh[jm] -= 1.0/h[Ny-1] * ( 0.5*v2*vb - vis4*vb/dy[Ny-1] );
		// mbc_2, j = 2
		v1 = 0.5 * (v[j1] + v[j2]);
		vb = vbc[j0];
		vis3 = nu[j0];
		rvh[j2] -= 1.0/h[2] * ( - 0.5*v1*vb - vis3*vb/dy[1] );

		// mbc_3, j = Ny-1
		v2 = 0.5 * (v[jn] + v[kmjn]);     w2 = w[jn];
		vb = 0.5 * (vbc[j1] + vbc[kmj1]); wb = wbc[j1];
		vis4 = 0.5 * (nu[jn] + nu[kmjn]);
		rwh[jm] -= 0.5/dy[Ny-1] * ( v2*wb + w2*vb - vis4*wb/h[Ny] - vis4*(v[jn]-v[kmjn])/dz );
		// mbc_3, j = 1
		v1 = 0.5 * (v[j1] + v[kmj1]);     w1 = w[j0];
		vb = 0.5 * (vbc[j0] + vbc[kmj0]); wb = wbc[j0];
		vis3 = 0.5 * (nu[j0] + nu[kmj0]);
		rwh[j1] -= 0.5/dy[1] * ( - v1*wb - w1*vb - vis3*wb/h[1] + vis3*(v[j1]-v[kmj1])/dz );
	}}
}






void IDM::getuh1(double *uh, double *u, double *v, double *w, double *nu)
/* compute deltaU^**, result returned by uh (the RHS ruh should be pre stored in uh ) */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, imjp, imkp, imjm, imkm;
	double u1, u2, v1, v2, w1, w2;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double *api = new double [Nx], *aci = new double [Nx], *ami = new double [Nx], *R1 = new double [Nx];
	double *apj = new double [Ny], *acj = new double [Ny], *amj = new double [Ny], *R2 = new double [Ny];
	double *apk = new double [Nz], *ack = new double [Nz], *amk = new double [Nz], *R3 = new double [Nz];
	Matrix matx1(Nx);
	Matrix maty1(Ny-1);
	Matrix matz1(Nz);

	# pragma omp for
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		for (j=1; j<Ny; j++) {
			idx = IDX(i,j,k);
			imjp = IDX(ima[i],j+1,k); imjm = IDX(ima[i],j-1,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

			v2 = 0.5 * ( v[jp] + v[imjp] );
			v1 = 0.5 * ( v[idx] + v[im] );
			vis3 = 0.25 * ( (nu[idx]+nu[im]) * dy[j-1] + (nu[jm]+nu[imjm]) * dy[j] ) / h[j];
			vis4 = 0.25 * ( (nu[idx]+nu[im]) * dy[j+1] + (nu[jp]+nu[imjp]) * dy[j] ) / h[j+1];

			apj[j] = ( 0.25/h[j+1]* v2                                   - 0.5/dy[j] * vis4/h[j+1]            ) * dt;
			acj[j] = ( 0.25/dy[j] *(v2*dy[j+1]/h[j+1] - v1*dy[j-1]/h[j]) + 0.5/dy[j] *(vis4/h[j+1]+vis3/h[j]) ) * dt + 1;
			amj[j] = (-0.25/h[j]  * v1                                   - 0.5/dy[j] * vis3/h[j]              ) * dt;

			R2 [j] = dt * uh[idx];
		}
		maty1.tdma( & amj[1], & acj[1], & apj[1], & R2[1] );// apj at j=Ny-1 and amj at j=1 are redundant in tdma, thus no need for explicit removal
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
		matx1.ctdma( ami, aci, api, R1 );
		for (i=0; i<Nx; i++) uh[IDX(i,j,k)] = R1[i];
	}}

	# pragma omp for
	for (j=1; j<Ny; j++) {
	for (i=0; i<Nx; i++) {
		for (k=0; k<Nz; k++) {
			idx = IDX(i,j,k);
			imkp = IDX(ima[i],j,kpa[k]); imkm = IDX(ima[i],j,kma[k]);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

			w2 = 0.5 * ( w[kp] + w[imkp] );
			w1 = 0.5 * ( w[idx] + w[im] );
			vis5 = 0.25 * ( nu[idx] + nu[im] + nu[km] + nu[imkm] );
			vis6 = 0.25 * ( nu[idx] + nu[im] + nu[kp] + nu[imkp] );

			apk[k] = ( 0.25/dz * w2     - 0.5/dz2 * vis6       ) * dt;
			ack[k] = ( 0.25/dz *(w2-w1) + 0.5/dz2 *(vis6+vis5) ) * dt + 1;
			amk[k] = (-0.25/dz * w1     - 0.5/dz2 * vis5       ) * dt;

			R3 [k] = uh[idx];
		}
		matz1.ctdma( amk, ack, apk, R3 );
		for (k=0; k<Nz; k++) uh[IDX(i,j,k)] = R3[k];
	}}

	delete [] api; delete [] aci; delete [] ami; delete [] R1;
	delete [] apj; delete [] acj; delete [] amj; delete [] R2;
	delete [] apk; delete [] ack; delete [] amk; delete [] R3;
}


void IDM::getuh2(double *uh, double *vh, double *u, double *v, double *w, double *nu)
/* compute deltaU^**, result returned by vh (the RHS rvh should be pre stored in uh ) */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, ipjm, jmkp, imjm, jmkm;
	double u1, u2, v1, v2, w1, w2, l21uh, m21uh;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double *api = new double [Nx], *aci = new double [Nx], *ami = new double [Nx], *R1 = new double [Nx];
	double *apj = new double [Ny], *acj = new double [Ny], *amj = new double [Ny], *R2 = new double [Ny];
	double *apk = new double [Nz], *ack = new double [Nz], *amk = new double [Nz], *R3 = new double [Nz];
	Matrix matx2(Nx);
	Matrix maty2(Ny-2);
	Matrix matz2(Nz);

	# pragma omp for
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		for (j=2; j<Ny; j++) {
			idx = IDX(i,j,k);
			ipjm = IDX(ipa[i],j-1,k); imjm = IDX(ima[i],j-1,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
		
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
		maty2.tdma( & amj[2], & acj[2], & apj[2], & R2[2] );
		for (j=2; j<Ny; j++) vh[IDX(i,j,k)] = R2[j];
	}}

	# pragma omp for
	for (j=2; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			idx = IDX(i,j,k);
			ipjm = IDX(ipa[i],j-1,k); imjm = IDX(ima[i],j-1,k);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

			u2 = ( u[ip]*dy[j-1] + u[ipjm]*dy[j] ) / (2.0*h[j]);
			u1 = ( u[idx]*dy[j-1] + u[jm]*dy[j] ) / (2.0*h[j]);
			vis1 = 0.25 * ( (nu[idx]+nu[im]) * dy[j-1] + (nu[jm]+nu[imjm]) * dy[j] ) / h[j];
			vis2 = 0.25 * ( (nu[idx]+nu[ip]) * dy[j-1] + (nu[jm]+nu[ipjm]) * dy[j] ) / h[j];

			api[i] = ( 0.25/dx * u2     - 0.5/dx2 * vis2       ) * dt;
			aci[i] = ( 0.25/dx *(u2-u1) + 0.5/dx2 *(vis2+vis1) ) * dt + 1;
			ami[i] = (-0.25/dx * u1     - 0.5/dx2 * vis1       ) * dt;

			R1 [i] = vh[idx];
		}
		matx2.ctdma( ami, aci, api, R1 );
		for (i=0; i<Nx; i++) vh[IDX(i,j,k)] = R1[i];
	}}

	# pragma omp for
	for (j=2; j<Ny; j++) {
	for (i=0; i<Nx; i++) {
		for (k=0; k<Nz; k++) {
			idx = IDX(i,j,k);
			jmkp = IDX(i,j-1,kpa[k]); jmkm = IDX(i,j-1,kma[k]);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

			w2 = ( w[kp]*dy[j-1] + w[jmkp]*dy[j] ) / (2.0*h[j]);
			w1 = ( w[idx]*dy[j-1] + w[jm]*dy[j] ) / (2.0*h[j]);
			vis5 = 0.25 * ( (nu[idx]+nu[km]) * dy[j-1] + (nu[jm]+nu[jmkm]) * dy[j] ) / h[j];
			vis6 = 0.25 * ( (nu[idx]+nu[kp]) * dy[j-1] + (nu[jm]+nu[jmkp]) * dy[j] ) / h[j];

			apk[k] = ( 0.25/dz * w2     - 0.5/dz2 * vis6       ) * dt;
			ack[k] = ( 0.25/dz *(w2-w1) + 0.5/dz2 *(vis6+vis5) ) * dt + 1;
			amk[k] = (-0.25/dz * w1     - 0.5/dz2 * vis5       ) * dt;

			R3 [k] = vh[idx];
		}
		matz2.ctdma( amk, ack, apk, R3 );
		for (k=0; k<Nz; k++) vh[IDX(i,j,k)] = R3[k];
	}}
	
	delete [] api; delete [] aci; delete [] ami; delete [] R1;
	delete [] apj; delete [] acj; delete [] amj; delete [] R2;
	delete [] apk; delete [] ack; delete [] amk; delete [] R3;
}


void IDM::getuh3(double *uh, double *vh, double *wh, double *u, double *v, double *w, double *nu)
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
	Matrix matx3(Nx);
	Matrix maty3(Ny-1);
	Matrix matz3(Nz);

	// ( I + dt M_33^2 )
	# pragma omp for
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		for (j=1; j<Ny; j++) {	jup = j==Ny-1 ? 0 : 1;	jum = j==1 ? 0 : 1;
			idx = IDX(i,j,k);
			ipkm = IDX(ipa[i],j,kma[k]); jpkm = IDX(i,j+1,kma[k]);
			imkm = IDX(ima[i],j,kma[k]); jmkm = IDX(i,j-1,kma[k]);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

			v2 = 0.5 * ( v[jp] + v[jpkm] );
			v1 = 0.5 * ( v[idx] + v[km] );
			vis3 = 0.25 * ( (nu[idx]+nu[km]) * dy[j-1] + (nu[jm]+nu[jmkm]) * dy[j] ) / h[j];
			vis4 = 0.25 * ( (nu[idx]+nu[km]) * dy[j+1] + (nu[jp]+nu[jpkm]) * dy[j] ) / h[j+1];

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
		maty3.tdma( & amj[1], & acj[1], & apj[1], & R2[1] );
		for (j=1; j<Ny; j++) wh[IDX(i,j,k)] = R2[j];
	}}

	// ( I + dt M_33^1 )
	# pragma omp for
	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			idx = IDX(i,j,k);
			ipkm = IDX(ipa[i],j,kma[k]); imkm = IDX(ima[i],j,kma[k]);
			ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
			im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

			u2 = 0.5 * ( u[ip] + u[ipkm] );
			u1 = 0.5 * ( u[idx] + u[km] );
			vis1 = 0.25 * ( nu[idx] + nu[im] + nu[km] + nu[imkm] );
			vis2 = 0.25 * ( nu[idx] + nu[ip] + nu[km] + nu[ipkm] );

			api[i] = ( 0.25/dx * u2     - 0.5/dx2 * vis2       ) * dt;
			aci[i] = ( 0.25/dx *(u2-u1) + 0.5/dx2 *(vis2+vis1) ) * dt + 1;
			ami[i] = (-0.25/dx * u1     - 0.5/dx2 * vis1       ) * dt;

			R1 [i] = wh[idx];
		}
		matx3.ctdma( ami, aci, api, R1 );
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
		matz3.ctdma( amk, ack, apk, R3 );
		for (k=0; k<Nz; k++) wh[IDX(i,j,k)] = R3[k];
	}}

	// update dvh
	# pragma omp for
	for (j=2; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		jmkp = IDX(i,j-1,kpa[k]); jmkm = IDX(i,j-1,kma[k]);
		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

		// m23wh
		v2 = 0.5 * ( v[idx] + v[kp] );
		v1 = 0.5 * ( v[idx] + v[km] );
		w2 = ( wh[kp]*dy[j-1] + wh[jmkp]*dy[j] ) / (2.0*h[j]);
		w1 = ( wh[idx]*dy[j-1] + wh[jm]*dy[j] ) / (2.0*h[j]);
		vis5 = 0.25 * ( nu[idx] + nu[jm] + nu[km] + nu[jmkm] );
		vis6 = 0.25 * ( nu[idx] + nu[jm] + nu[kp] + nu[jmkp] );

		l23wh = ( vis6 * (wh[kp]-wh[jmkp]) - vis5 * (wh[idx]-wh[jm]) ) / (2.0*h[j]*dz);
		m23wh = (v2*w2 - v1*w1) / (2.0*dz) - l23wh;

		vh[idx] -= dt * m23wh;
	}}}
	
	// update duh
	# pragma omp for
	for (j=1; j<Ny; j++) {	jup = j==Ny-1 ? 0 : 1;	jum = j==1 ? 0 : 1;
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		imjp = IDX(ima[i],j+1,k); imkp = IDX(ima[i],j,kpa[k]);
		imjm = IDX(ima[i],j-1,k); imkm = IDX(ima[i],j,kma[k]);
		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);

		// m12vh
		u2 = ( u[idx]*dy[j+1] + u[jp]*dy[j] ) / (2.0*h[j+1]);
		u1 = ( u[idx]*dy[j-1] + u[jm]*dy[j] ) / (2.0*h[j]);
		v2 = 0.5 * ( vh[jp] + vh[imjp] );
		v1 = 0.5 * ( vh[idx] + vh[im] );
		vis3 = 0.25 * ( (nu[idx]+nu[im]) * dy[j-1] + (nu[jm]+nu[imjm]) * dy[j] ) / h[j];
		vis4 = 0.25 * ( (nu[idx]+nu[im]) * dy[j+1] + (nu[jp]+nu[imjp]) * dy[j] ) / h[j+1];

		u2 *= jup;	vis4 *= jup;
		u1 *= jum;	vis3 *= jum;
		l12vh = ( vis4 * (vh[jp]-vh[imjp]) - vis3 * (vh[idx]-vh[im]) ) / (2.0*dx*dy[j]);
		m12vh = (u2*v2 - u1*v1) / (2.0*dy[j]) - l12vh;

		// m13wh
		u2 = 0.5 * ( u[idx] + u[kp] );
		u1 = 0.5 * ( u[idx] + u[km] );
		w2 = 0.5 * ( wh[kp] + wh[imkp] );
		w1 = 0.5 * ( wh[idx] + wh[im] );
		vis5 = 0.25 * ( nu[idx] + nu[im] + nu[km] + nu[imkm] );
		vis6 = 0.25 * ( nu[idx] + nu[im] + nu[kp] + nu[imkp] );

		l13wh = ( vis6 * (wh[kp]-wh[imkp]) - vis5 * (wh[idx]-wh[im]) ) / (2.0*dx*dz);
		m13wh = (u2*w2 - u1*w1) / (2.0*dz) - l13wh;

		uh[idx] -= dt * ( m12vh + m13wh );
	}}}

	// update intermediate velocity field
	# pragma omp for
	for (j=1; j<Ny; j++) {	jum = j==1 ? 0 : 1;
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		uh[idx] += u[idx];
		vh[idx] += v[idx] * jum; // j=1 is redundant
		wh[idx] += w[idx];
	}}}

	delete [] api; delete [] aci; delete [] ami; delete [] R1;
	delete [] apj; delete [] acj; delete [] amj; delete [] R2;
	delete [] apk; delete [] ack; delete [] amk; delete [] R3;
}





void IDM::rhsdp(double *rdp, double *uh, double *vh, double *wh, double *vbc)
{
	int i, j, k, idx, ip, jp, kp, jbp, jbm, jup, jum;
	for (j=1; j<Ny; j++) {	jup = j==Ny-1 ? 0 : 1;	jum = j==1 ? 0 : 1;
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);	ip = IDX(ipa[i],j,k);
		jbp = IDX(i,1,k);	jp = IDX(i,j+1,k);
		jbm = IDX(i,0,k);	kp = IDX(i,j,kpa[k]);
		// ( Du^h - cbc ) / dt
		rdp[idx] = 1.0/dt * (
			( uh[ip] - uh[idx] ) / dx
		+	( wh[kp] - wh[idx] ) / dz
		+	( vh[jp] * jup + vbc[jbp] * (1-jup)
			- vh[idx]* jum - vbc[jbm] * (1-jum) ) / dy[j]	);
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

		// set reference pressure P(kx=0,kz=0,j=1) to be 0
		if (k==0 && i==0) {
			cpj[1] = 0;
			ccj[1] = 1;
			cmj[1] = 0;
			cfj1[1] = - Nxz * refp;	// fdp(kx=0,kz=0,j=1) = -Nxz * fp(kx=0,kz=0,j=1) ==> <dp(j=1)> = - <p(j=1)>
			cfj2[1] = 0;
		}

		// if (k==0 && i==0) ccj[1] += 1e20;

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






void IDM::update(
	double *u,	double *v,	double *w,	double *p,
	double *uh,	double *vh,	double *wh,	double *dp)
/* project from U^* to U^n+1 using DP */
{
	int i, j, k, idx, im, jm, km;
	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = IDX(i,j,k);
		im = IDX(ima[i],j,k);
		jm = IDX(i,j-1,k);
		km = IDX(i,j,kma[k]);
		// store time derivative in intermediate variables
		uh[idx] = (uh[idx] - u[idx]) / dt - (dp[idx] - dp[im]) / dx;
		if (j>1) vh[idx] = (vh[idx] - v[idx]) / dt - (dp[idx] - dp[jm]) / h[j];
		wh[idx] = (wh[idx] - w[idx]) / dt - (dp[idx] - dp[km]) / dz;
		// update fields
		u[idx] += uh[idx] * dt;
		if (j>1) v[idx] += vh[idx] * dt;
		w[idx] += wh[idx] * dt;
		p[idx] += dp[idx];
	}}}

	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) { dp[IDX(i,j,k)] /= dt; }}}
}


void IDM::meanpg(double *mpg, Vctr &U)
/* solve the increment of mean pressure gradient at n+1/2 step, given mass flow rate at n+1 step */
{
	int i, j, k;
	Scla &U1 = U.com1, &U3 = U.com3;
	// mean pressure gradient increment is solved by fixing streamwise flow rate 2.0 and spanwise flow rate 0
	double dmpg1 = (U1.bulkMeanU() - 1.0) / dt;
	double dmpg3 =  U3.bulkMeanU()        / dt;
	// update the mean pressure gradient
	mpg[0] += dmpg1;
	mpg[2] += dmpg3;
	// complement the mean pressure gradient increment that was not included in the velocity update step
	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		U1.id(i,j,k) -= dt * dmpg1;
		U3.id(i,j,k) -= dt * dmpg3;
	}}}
}












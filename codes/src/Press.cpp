# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>

# include "Press.h"
# include "Matrix.h"



void Press::rhs(Scla &R, Scla &NU, Vctr &U, Vctr &V, double dt)
/* compute the RHS of Poisson equation at step n+1/2. V is the velocity field at step n+1, and U is that at step n */
{
	Mesh &ms = R.meshGet();
	int Nx = ms.Nx, Ny = ms.Ny, Nz = ms.Nz;
	double dx = ms.dx, dz = ms.dz, dx2 = ms.dx2, dz2 = ms.dz2;
	double *dy = ms.dy, *h = ms.h;
	double *hm = ms.hm, *hc = ms.hc, *hp = ms.hp;
	double *dym= ms.dym,*dyc= ms.dyc,*dyp= ms.dyp;
	int *kpa = ms.kpa, *kma = ms.kma;
	int *ipa = ms.ipa, *ima = ms.ima;

	int i, j, k, idx, ip, im, jp, jm, kp, km, jup, jum;
	int imjp, imkp, imjm, imkm, ipjm, jmkp, jmkm, ipkm, jpkm;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double u1, u2, v1, v2, w1, w2;
	double api, aci, ami, apj, acj, amj, apk, ack, amk;
	double l11un, l12vn, l13wn, m11un, m12vn, m13wn;
	double l21un, l22vn, l23wn, m21un, m22vn, m23wn;
	double l31un, l32vn, l33wn, m31un, m32vn, m33wn;

	Vctr RUH(ms);
	double *ruh, *rvh, *rwh, *u, *v, *w, *nu;

	nu = NU.bulkGet();
	ruh = RUH.com1.bulkGet();
	rvh = RUH.com2.bulkGet();
	rwh = RUH.com3.bulkGet();

	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {

		idx = ms.IDX(i,j,k);
		ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
		im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);
		imjp = ms.IDX(ima[i],j+1,k); imkp = ms.IDX(ima[i],j,kpa[k]);
		imjm = ms.IDX(ima[i],j-1,k); imkm = ms.IDX(ima[i],j,kma[k]);
		ipjm = ms.IDX(ipa[i],j-1,k); ipkm = ms.IDX(ipa[i],j,kma[k]);
		jmkm = ms.IDX(i,j-1,kma[k]); jmkp = ms.IDX(i,j-1,kpa[k]);
		jpkm = ms.IDX(i,j+1,kma[k]);
	
		//// U component
		// interpolate viscosity and velocity at step n to the position needed
		vis1 = nu[im];
		vis2 = nu[idx];
		vis3 = 0.25 * ( (nu[idx]+nu[im]) * dy[j-1] + (nu[jm]+nu[imjm]) * dy[j] ) / h[j];
		vis4 = 0.25 * ( (nu[idx]+nu[im]) * dy[j+1] + (nu[jp]+nu[imjp]) * dy[j] ) / h[j+1];
		vis5 = 0.25 * ( nu[idx] + nu[im] + nu[km] + nu[imkm] );
		vis6 = 0.25 * ( nu[idx] + nu[im] + nu[kp] + nu[imkp] );

		u = U.com1.bulkGet(); v = U.com2.bulkGet(); w = U.com3.bulkGet();
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

		u = V.com1.bulkGet();
		m11un =	api*u[ip] + aci*u[idx] + ami*u[im]
			+	apj*u[jp] + acj*u[idx] + amj*u[jm]
			+	apk*u[kp] + ack*u[idx] + amk*u[km];

		// m12vn
		u = U.com1.bulkGet();
		u2 = ( u[idx]*dy[j+1] + u[jp]*dy[j] ) / (2.0*h[j+1]);
		u1 = ( u[idx]*dy[j-1] + u[jm]*dy[j] ) / (2.0*h[j]);

		v = V.com2.bulkGet();
		v2 = 0.5 * (v[jp] + v[imjp]);
		v1 = 0.5 * (v[idx]+ v[im]);

		m12vn = (u2*v2 - u1*v1) / (2.0*dy[j]) - ( vis4 * (v[jp]-v[imjp]) - vis3 * (v[idx]-v[im]) ) / (2.0*dx*dy[j]);

		// m13wn
		u = U.com1.bulkGet();
		u2 = 0.5 * ( u[idx] + u[kp] );
		u1 = 0.5 * ( u[idx] + u[km] );

		w = V.com3.bulkGet();
		w2 = 0.5 * (w[kp] + w[imkp]);
		w1 = 0.5 * (w[idx]+ w[im]);

		m13wn = (u2*w2 - u1*w1) / (2.0*dz) - ( vis6 * (w[kp]-w[imkp]) - vis5 * (w[idx]-w[im]) ) / (2.0*dx*dz);

		// R_1 without boundary modification
		ruh[idx] = (l11un + l12vn + l13wn) - (m11un + m12vn + m13wn);


		//// W component
		// interpolate viscosity and velocity at step n to the position needed
		vis1 = 0.25 * ( nu[idx] + nu[im] + nu[km] + nu[imkm] );
		vis2 = 0.25 * ( nu[idx] + nu[ip] + nu[km] + nu[ipkm] );
		vis3 = 0.25 * ( (nu[idx]+nu[km]) * dy[j-1] + (nu[jm]+nu[jmkm]) * dy[j] ) / h[j];
		vis4 = 0.25 * ( (nu[idx]+nu[km]) * dy[j+1] + (nu[jp]+nu[jpkm]) * dy[j] ) / h[j+1];
		vis5 = nu[km];
		vis6 = nu[idx];

		u = U.com1.bulkGet(); v = U.com2.bulkGet(); w = U.com3.bulkGet();
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

		w = V.com3.bulkGet();
		m33wn =	api*w[ip] + aci*w[idx] + ami*w[im]
			+	apj*w[jp] + acj*w[idx] + amj*w[jm]
			+	apk*w[kp] + ack*w[idx] + amk*w[km];

		// m32vn
		w = U.com3.bulkGet();
		w2 = ( w[idx]*dy[j+1] + w[jp]*dy[j] ) / (2.0*h[j+1]);
		w1 = ( w[idx]*dy[j-1] + w[jm]*dy[j] ) / (2.0*h[j]);

		v = V.com2.bulkGet();
		v2 = 0.5 * (v[jp] + v[jpkm]);
		v1 = 0.5 * (v[idx]+ v[km]);

		m32vn = (w2*v2 - w1*v1) / (2.0*dy[j]) - ( vis4 * (v[jp]-v[jpkm]) - vis3 * (v[idx]-v[km]) ) / (2.0*dz*dy[j]);

		// m31un
		w = U.com3.bulkGet();
		w2 = 0.5 * ( w[idx] + w[ip] );
		w1 = 0.5 * ( w[idx] + w[im] );

		u = V.com1.bulkGet();
		u2 = 0.5 * (u[ip] + u[ipkm]);
		u1 = 0.5 * (u[idx]+ u[km]);

		m31un = (w2*u2 - w1*u1) / (2.0*dx) - ( vis2 * (u[ip]-u[ipkm]) - vis1 * (u[idx]-u[km]) ) / (2.0*dz*dx);

		// R_3 without boundary modification
		rwh[idx] = (l31un + l32vn + l33wn) - (m31un + m32vn + m33wn);

		if (j == 1) continue;

		//// V component
		// interpolate viscosity and velocity at step n to the position needed
		vis1 = 0.25 * ( (nu[idx]+nu[im]) * dy[j-1] + (nu[jm]+nu[imjm]) * dy[j] ) / h[j];
		vis2 = 0.25 * ( (nu[idx]+nu[ip]) * dy[j-1] + (nu[jm]+nu[ipjm]) * dy[j] ) / h[j];
		vis3 = nu[jm];
		vis4 = nu[idx];
		vis5 = 0.25 * ( (nu[idx]+nu[km]) * dy[j-1] + (nu[jm]+nu[jmkm]) * dy[j] ) / h[j];
		vis6 = 0.25 * ( (nu[idx]+nu[kp]) * dy[j-1] + (nu[jm]+nu[jmkp]) * dy[j] ) / h[j];

		u = U.com1.bulkGet(); v = U.com2.bulkGet(); w = U.com3.bulkGet();
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

		v = V.com2.bulkGet();
		m22vn =	api*v[ip] + aci*v[idx] + ami*v[im]
			+	apj*v[jp] + acj*v[idx] + amj*v[jm]
			+	apk*v[kp] + ack*v[idx] + amk*v[km];

		// m21un
		v = U.com2.bulkGet();
		v2 = 0.5 * ( v[idx] + v[ip] );
		v1 = 0.5 * ( v[idx] + v[im] );

		u = V.com1.bulkGet();
		u2 = ( u[ip]*dy[j-1] + u[ipjm]*dy[j] ) / (2.0*h[j]);
		u1 = ( u[idx]*dy[j-1] + u[jm]*dy[j] ) / (2.0*h[j]);

		m21un = (v2*u2 - v1*u1) / (2.0*dx) - ( vis2 * (u[ip]-u[ipjm]) - vis1 * (u[idx]-u[jm]) ) / (2.0*h[j]*dx);

		// m23wn
		v = U.com2.bulkGet();
		v2 = 0.5 * ( v[idx] + v[kp] );
		v1 = 0.5 * ( v[idx] + v[km] );

		w = V.com3.bulkGet();
		w2 = ( w[kp]*dy[j-1] + w[jmkp]*dy[j] ) / (2.0*h[j]);
		w1 = ( w[idx]*dy[j-1] + w[jm]*dy[j] ) / (2.0*h[j]);

		m23wn = (v2*w2 - v1*w1) / (2.0*dz) - ( vis6 * (w[kp]-w[jmkp]) - vis5 * (w[idx]-w[jm]) ) / (2.0*h[j]*dz);

		// R_2 without boundary modification
		rvh[idx] = (l21un + l22vn + l23wn) - (m21un + m22vn + m23wn);

	}}}


	u = U.com2.bulkGet();
	v = V.com2.bulkGet();

	for (j=1; j<Ny; j++) { jup = (j!=Ny-1); jum = (j!=1);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx= ms.IDX(i,j,k);
		ip = ms.IDX(ipa[i],j,k);
		jp = ms.IDX(i,j+1,k);
		kp = ms.IDX(i,j,kpa[k]);
		R.id(i,j,k) =
			( ruh[ip] - ruh[idx] ) / dx
		+	( rwh[kp] - rwh[idx] ) / dz
		+	( rvh[jp] * jup - rvh[idx] * jum ) / dy[j]
		+	( v[jp]*(1-jup) - v[idx]*(1-jum) ) / dy[j] / dt
		-	( u[jp]*(1-jup) - u[idx]*(1-jum) ) / dy[j] / dt;
	}}}
}

void Press::poisson(Scla &P)
/* solve the Poisson equation, with RHS pre-stored in P */
{
	Mesh &ms = P.meshGet();
	double *ak1 = ms.ak1, *ak3 = ms.ak3;
	double *ppj = ms.ppj, *pcj = ms.pcj, *pmj = ms.pmj;
	int Nxc = ms.Nx/2+1, Nxr = 2*Nxc, Nxzr = ms.Nz*Nxr, Ny = ms.Ny, Nz = ms.Nz;

	int i, j, k;
	double *cpj = new double [Ny], *ccj = new double [Ny], *cmj = new double [Ny];
	double *cfj1= new double [Ny], *cfj2= new double [Ny];
	Matrix matyp(Ny-1);
	
	P.fft();

	for (k=0; k<Nz; k++) {
	for (i=0; i<Nxc; i++) {
		for (j=1; j<Ny; j++) {
			cpj[j] = ppj[j];
			ccj[j] = pcj[j] - ak3[k] - ak1[i];
			cmj[j] = pmj[j];
			cfj1[j]= P.idf(2*i,  j,k);	// real part
			cfj2[j]= P.idf(2*i+1,j,k);	// imaginary part
		}
		// set reference pressure P(kx=0,kz=0,j=1) to be 0
		if (k==0 && i==0) {
			cpj[1] = 0;
			ccj[1] = 1;
			cmj[1] = 0;
			cfj1[1]= 0;
			cfj2[1]= 0;
		}

		matyp.tdma( & cmj[1], & ccj[1], & cpj[1], & cfj1[1] );
		matyp.tdma( & cmj[1], & ccj[1], & cpj[1], & cfj2[1] );

		for (j=1; j<Ny; j++) {
			P.idf(2*i,  j,k) = cfj1[j];
			P.idf(2*i+1,j,k) = cfj2[j];
		}
	}}

	P.ifft();

	delete [] cpj; delete [] ccj; delete [] cmj;
	delete [] cfj1;delete [] cfj2;
}



// #define DEBUG
#ifdef DEBUG

void main()
{
	Mesh ms(256,129,256,12.5664,2,6.2832,0,"/run/media/root/DATA/whn/channel_IDM/data_DNS180/statdata/");
	Vctr U(ms), UT(ms);
	Scla NU(ms), P0(ms), P(ms, true);
	Press pres;
	double dt = 5e-3;

	NU.bulkSet(1./2850);
	P0.fileIO("/run/media/root/DATA/whn/channel_IDM/data_DNS180/fielddata/", "P00025000", 'r');
	U.com1.fileIO("/run/media/root/DATA/whn/channel_IDM/data_DNS180/fielddata/", "U00025000", 'r');
	U.com2.fileIO("/run/media/root/DATA/whn/channel_IDM/data_DNS180/fielddata/", "V00025000", 'r');
	U.com3.fileIO("/run/media/root/DATA/whn/channel_IDM/data_DNS180/fielddata/", "W00025000", 'r');
	UT.com1.fileIO("/run/media/root/DATA/whn/channel_IDM/data_DNS180/fielddata/", "UT00025000", 'r');
	UT.com2.fileIO("/run/media/root/DATA/whn/channel_IDM/data_DNS180/fielddata/", "VT00025000", 'r');
	UT.com3.fileIO("/run/media/root/DATA/whn/channel_IDM/data_DNS180/fielddata/", "WT00025000", 'r');

	pres.U2P(P, NU, U, UT, dt);

	P.debug_AsciiOutput("", "P", 1, 129);
	P0.debug_AsciiOutput("", "P0", 1, 129);
}



/***** example of Makefile *****/

// ##### path names #####
// # header files, source files, middle files and data files
// # relative to Makefile

// IDIR = include
// SDIR = src
// ODIR = obj
// DATADIR = ../data


// ##### compile & link options #####

// OUTPUT = press_test

// CC = icc
// LIBS = -lfftw3 -lm
// CFLAGS = -I $(IDIR) $(LIBS)


// ##### file lists #####

// _DEPS = Basic.h Matrix.h Press.h
// DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

// _OBJS = Bulk.o Mesh.o Scla.o Vctr.o Matrix.o Press.o
// OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))


// #### targets & commands #####
// # note: $@ target, $^ all dependents, $< the first dependent

// $(OUTPUT): $(OBJS)
// 	$(CC) -o $@ $^ $(CFLAGS)


// $(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
// 	$(CC) -c -o $@ $< $(CFLAGS)

// .PHONY: remake
// remake:
// 	make clean
// 	make
// 	mv $(OUTPUT) $(DATADIR)/

// .PHONY: clean
// clean:
// 	rm -f $(OUTPUT) *~ $(ODIR)/*.o $(INCDIR)/*~





/***** example of check file *****/

// #!/root/Software/anaconda3/bin/python3
// import numpy as np
// import matplotlib
// matplotlib.use("Agg")
// import matplotlib.pyplot as plt


// P = np.loadtxt("P.txt").reshape([128, 256, 256])
// P0 = np.loadtxt("P0.txt").reshape([128, 256, 256])


// print(np.mean(abs(P-P0)))


// fig, axs = plt.subplots(2,1)
// axs[0].contour(P[0])
// axs[1].contour(P0[0])
// fig.savefig("test1.png")


// Pm = np.mean(P, axis=(-1,-2))
// P0m = np.mean(P0, axis=(-1,-2))
// P -= Pm.reshape([128,1,1])
// P0 -= P0m.reshape([128,1,1])

// fig, axs = plt.subplots(1,2)
// axs[0].plot(Pm)
// axs[0].plot(P0m)
// axs[1].plot(np.mean(P**2, axis=(-1,-2)))
// axs[1].plot(np.mean(P0**2, axis=(-1,-2)))
// fig.savefig("test2.png")





/***** example of dynamic library Makefile *****/

// ##### path names #####
// # header files, source files, middle files and data files
// # relative to Makefile

// IDIR = include
// SDIR = src
// ODIR = obj
// DATADIR = ../data


// ##### compile & link options #####

// OUTPUT = libpress.so

// CC = icc
// LIBS = -lfftw3 -lm
// CFLAGS = -I $(IDIR) $(LIBS)


// ##### file lists #####

// _DEPS = Basic.h Matrix.h Press.h
// DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

// _OBJS = Bulk.o Mesh.o Scla.o Vctr.o Matrix.o Press.o
// OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))


// #### targets & commands #####
// # note: $@ target, $^ all dependents, $< the first dependent

// $(OUTPUT): $(OBJS)
// 	$(CC) -shared -o $@ $^ $(CFLAGS)


// $(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
// 	$(CC) -fPIC -c -o $@ $< $(CFLAGS)

// .PHONY: remake
// remake:
// 	make clean
// 	make
// 	mv $(OUTPUT) $(DATADIR)/

// .PHONY: clean
// clean:
// 	rm -f $(OUTPUT) *~ $(ODIR)/*.o $(INCDIR)/*~


#endif

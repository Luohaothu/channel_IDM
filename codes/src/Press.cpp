# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>

# include "Press.h"
# include "Matrix.h"
# include "IDM.h"

using namespace std;




void Press::rhs(Feld &R, const Vctr &U, const Vctr &UT, const Feld &VIS, double dt)
{
	Mesh mesh(R.meshGet());
	Mesh bmesh(mesh.Nx,2-1,mesh.Nz, mesh.Lx,mesh.Ly,mesh.Lz);

	bmesh.y[0] = 0; bmesh.yc[0] = mesh.y[1];
	bmesh.y[1] = 0; bmesh.yc[1] = mesh.y[mesh.Ny];

	Feld FLD(mesh), BC(bmesh); IDM idm(mesh);

	// set FLD to u^n and p=0
	((FLD.V[1] = UT[1]) *= -dt) += U[1];
	((FLD.V[2] = UT[2]) *= -dt) += U[2];
	((FLD.V[3] = UT[3]) *= -dt) += U[3];
	FLD.S = 0.;
	// get velocity BC from u^n+1
	BC.V[1].lyrSet(U[1][0], 0); BC.V[1].lyrSet(U[1][mesh.Ny], 1);
	BC.V[2].lyrSet(U[2][1], 0); BC.V[2].lyrSet(U[2][mesh.Ny], 1);
	BC.V[3].lyrSet(U[3][0], 0); BC.V[3].lyrSet(U[3][mesh.Ny], 1);
	// solve for RHS of Gp^n+.5
	R.reset();
	idm.urhs(R.V, FLD, VIS, BC, R.V); // R.V act as FB and is pre-set to 0
	idm.muh1(FLD.S.blkGet(), UT, FLD.V, VIS); R.V[1] -= ((FLD.S *= dt) += UT[1]);
	idm.muh2(FLD.S.blkGet(), UT, FLD.V, VIS); R.V[2] -= ((FLD.S *= dt) += UT[2]);
	idm.muh3(FLD.S.blkGet(), UT, FLD.V, VIS); R.V[3] -= ((FLD.S *= dt) += UT[3]);
	// pre-eliminate boundary before operated by operator D
	R.V[1].lyrSet(0., 0); R.V[1].lyrSet(0., mesh.Ny);
	R.V[2].lyrSet(0., 1); R.V[2].lyrSet(0., mesh.Ny);
	R.V[3].lyrSet(0., 0); R.V[3].lyrSet(0., mesh.Ny);
	// solve for RHS of DGp^n+.5
	for (int j=1; j<mesh.Ny; j++) {
	for (int k=0; k<mesh.Nz; k++) {
	for (int i=0; i<mesh.Nx; i++) {
		R.S.id(i,j,k) = R.V.divergence(i,j,k);
	}}}

	bmesh.freeall();
}

void Press::poisson(Scla &P)
/* solve the Poisson equation, with RHS pre-stored in P */
{
	const Mesh &ms = P.meshGet();
	int i, j, k, Nxc = ms.Nx/2+1, Ny = ms.Ny, Nz = ms.Nz;
	double *cpj = new double [Ny];
	double *ccj = new double [Ny];
	double *cmj = new double [Ny];
	double *cfj1= new double [Ny], *cfj2= new double [Ny];
	Matrix matyp(Ny-1);
	
	P.fft();

	for (k=0; k<Nz; k++) {
	for (i=0; i<Nxc; i++) {
		for (j=1; j<Ny; j++) {
			cpj[j] = ms.ppj[j];
			ccj[j] = ms.pcj[j] - ms.ak3[k] - ms.ak1[i];
			cmj[j] = ms.pmj[j];
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

int main()
{
	Mesh ms(32,49,32, 6.2832,2,3.1416, 0,"/Users/whn/Desktop/test/data1/statdata/");
	Feld FLD(ms), VIS(ms); Vctr UT(ms);
	double dt = 5e-3;

	VIS.reset(1./2850);

	FLD.S.fileIO("/Users/whn/Desktop/test/data1/fielddata/", "P00005000", 'r');
	FLD.S.debug_AsciiOutput("", "P0", 0, 50);

	FLD.V[1].fileIO("/Users/whn/Desktop/test/data1/fielddata/", "U00005000", 'r');
	FLD.V[2].fileIO("/Users/whn/Desktop/test/data1/fielddata/", "V00005000", 'r');
	FLD.V[3].fileIO("/Users/whn/Desktop/test/data1/fielddata/", "W00005000", 'r');
	UT[1].fileIO("/Users/whn/Desktop/test/data1/fielddata/", "UT00005000", 'r');
	UT[2].fileIO("/Users/whn/Desktop/test/data1/fielddata/", "VT00005000", 'r');
	UT[3].fileIO("/Users/whn/Desktop/test/data1/fielddata/", "WT00005000", 'r');

	Press::u2p(FLD, UT, VIS, dt);
	FLD.S.debug_AsciiOutput("", "P", 0, 50);

	cout << "P written." << endl;
	return 0;
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

// # CC = icc
// # LIBS = -lfftw3 -lm -qopenmp
// CC = g++-9
// LIBS = -lfftw3 -lm -fopenmp
// CFLAGS = -I $(IDIR) $(LIBS)


// ##### file lists #####

// _DEPS = Basic.h Interp.h Matrix.h IDM.h Press.h
// DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

// _OBJS = Mesh.o Bulk.o Scla.o Vctr.o Feld.o Interp.o Matrix.o IDM.o Press.o
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


// P = np.loadtxt("P.txt").reshape([32, 50, 32])
// P0 = np.loadtxt("P0.txt").reshape([32, 50, 32])


// print(np.mean(abs(P-P0)))


// fig, axs = plt.subplots(2,1)
// axs[0].contour(P[0])
// axs[1].contour(P0[0])
// fig.savefig("test1.png")


// Pm = np.mean(P, axis=(-1,-2))
// P0m = np.mean(P0, axis=(-1,-2))
// P -= Pm.reshape([32,1,1])
// P0 -= P0m.reshape([32,1,1])

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





# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>

# include "Interp.h"

using namespace std;


Interp::Interp(const Scla &src, Scla &dst):
src(src), ms0(src.meshGet()),
dst(dst), ms1(dst.meshGet()),
Nx0(ms0.Nx), Nx1(ms1.Nx),
Ny0(ms0.Ny), Ny1(ms1.Ny),
Nz0(ms0.Nz), Nz1(ms1.Nz),
Nxz0(ms0.Nxz),Nxz1(ms1.Nxz),
Lx0(ms0.Lx), Lx1(ms1.Lx),
Ly0(ms0.Ly), Ly1(ms1.Ly),
Lz0(ms0.Lz), Lz1(ms1.Lz),
dx0(Lx0/Nx0),dx1(Lx1/Nx1),
dz0(Lz0/Nz0),dz1(Lz1/Nz1),
rscl0(1.0*Nxz1/Nxz0),
rscl1(Lx0*Lz0/Lx1/Lz1)
{
	if (dx0>=Lx1 || dx1>=Lx0 || dz0>=Lz1 || dz1>=Lz0)
	{ cout << "Invalid mesh for interpolation !" << endl; exit(0); }
	if (Ly1 > Ly0)
	{ cout << "No extrapolation allowed for y mesh !" << endl; exit(0); }

	i0l = new int [Nx1];
	k0l = new int [Nz1];

	i0f = new int [Nx1];
	k0f = new int [Nz1];
	nif = new int [Nx1];
	nkf = new int [Nz1];

	j0u = new int [Ny1+1];
	j0v = new int [Ny1+1];
	j0x = new int [Ny1+1];

	int i0, j0, k0, i1, j1, k1;
	double x0, y0, z0, x1, y1, z1;

	// for every X1, find the largest X0 that is <= X1
	for (i1=0; i1<Nx1; i1++) {
		x1 = i1 * dx1;
		if (i1 == 0) x0 = (i0=1) * dx0; // every search begins from the last position
		while (x0 <= x1) x0 = (++i0) * dx0;
		i0l[i1] = i0 - 1;
	}
	// for every Z1, find the largest Z0 that is <= Z1
	for (k1=0; k1<Nz1; k1++) {
		z1 = k1 * dz1;
		if (k1 == 0) z0 = (k0=1) * dz0;
		while (z0 <= z1) z0 = (++k0) * dz0;
		k0l[k1] = k0 - 1;
	}

	// for every X1, find the largest X0 that is <= X1-0.5dx1, and the smallest X0 that is >= X1+0.5dx1
	for (i1=0; i1<Nx1; i1++) {
		x1 = (i1-0.5) * dx1;
		if (i1 == 0) x0 = (i0=1-Nx0) * dx0;
		while (x0 <= x1) x0 = (++i0) * dx0;
		i0f[i1] = i0 - 1;

		for (nif[i1]=1; (nif[i1]+i0f[i1])*dx0 < x1+dx1; nif[i1]++) {}
	}
	// for every Z1, find the largest Z0 that is <= Z1-0.5dz1, and the smallest Z0 that is >= Z1+0.5dz1
	for (k1=0; k1<Nz1; k1++) {
		z1 = (k1-0.5) * dz1;
		if (k1 == 0) z0 = (k0=1-Nz0) * dz0;
		while (z0 <= z1) z0 = (++k0) * dz0;
		k0f[k1] = k0 - 1;

		for (nkf[k1]=1; (nkf[k1]+k0f[k1])*dz0 < z1+dz1; nkf[k1]++) {}
	}

	// for every Y1, find the largest Y0 that is <= Y1
	// if Y1 is out of Y0 range, the value will be linearly extrapolated
	// yc -> yc
	for (j1=0; j1<=Ny1; j1++) {
		y1 = ms1.yc[j1];
		if (j1 == 0) y0 = ms0.yc[j0=1];
		while (y0 <= y1 && j0 < Ny0) y0 = ms0.yc[++j0];
		j0u[j1] = j0 - 1;
	}
	// y -> y
	j0v[0] = 0; // note: Y1[0] do not participate in interpolation
	for (j1=1; j1<=Ny1; j1++) {
		y1 = ms1.y[j1];
		if (j1 == 1) y0 = ms0.y[j0=2];
		while (y0 <= y1 && j0 < Ny0) y0 = ms0.y[++j0];
		j0v[j1] = j0 - 1;
	}
	// y -> yc
	for (j1=0; j1<=Ny1; j1++) {
		y1 = ms1.yc[j1];
		if (j1 == 0) y0 = ms0.y[j0=2];
		while (y0 <= y1 && j0 < Ny0) y0 = ms0.y[++j0];
		j0x[j1] = j0 - 1;
	}
}

Interp::~Interp()
{
	delete [] i0l; delete [] k0l;
	delete [] i0f; delete [] k0f;
	delete [] nif; delete [] nkf;
	delete [] j0u; delete [] j0v;
}



void Interp::layerPrdLin(int j0, int j1)
/* planar linear interpolation in X & Z direction */
/* from j0 layer of src to j1 layer of dst, periodically extrapolated if needed */
{
	int i, k, im0, km0, ip0, kp0;
	double x, z, x1, z1, x2, z2;

	// bilinear interpolation
	for (k=0; k<Nz1; k++) {

		z = k * dz1;
		z1 = k0l[k] * dz0; z2 = (k0l[k]+1) * dz0;
		km0= k0l[k] % Nz0; kp0= (k0l[k]+1) % Nz0;

		for (i=0; i<Nx1; i++) {

			x = i * dx1;
			x1 = i0l[i] * dx0; x2 = (i0l[i]+1) * dx0;
			im0= i0l[i] % Nx0; ip0= (i0l[i]+1) % Nx0;

			dst.id(i,j1,k) = 1.0 / (x2-x1) / (z2-z1) * (
				(x2-x) * (z2-z) * src.id(im0,j0,km0) + \
				(x-x1) * (z2-z) * src.id(ip0,j0,km0) + \
				(x2-x) * (z-z1) * src.id(im0,j0,kp0) + \
				(x-x1) * (z-z1) * src.id(ip0,j0,kp0)	);
	}}
}



void Interp::layerPrdFlt(int j0, int j1)
/* planar filter in X & Z direction */
/* from j0 layer of src to j1 layer of dst, periodically extrapolated if needed */
/* filter is defined as the integration of src in a cell determined by the dst grid */
{
	int i, k, i0, k0, ii, kk;
	double x, z, x0, z0;
	double delta, weight, wx, wz, integr;

	for (k=0; k<Nz1; k++) { z = k * dz1;
	for (i=0; i<Nx1; i++) { x = i * dx1;

		weight = 0;
		integr = 0;

		// trapezoidal integration in two directions respectively
		for (kk=0; kk<=nkf[k]; kk++) {

			z0 = (k0f[k]+kk) * dz0;
			k0 = (k0f[k]+kk+Nx0) % Nz0;
	
			if      ((delta = z-dz1/2.-z0) > dz0) wz = 0; // not likely
			else if (delta > 0)                   wz = .5 * (1-delta/dz0) * (dz0-delta);
			else if (delta > -dz0)                wz = .5 * (-delta/dz0) * (2*dz0+delta) + dz0/2.;
			else if ((delta = z+dz1/2.-z0) <-dz0) wz = 0; // not likely
			else if (delta < 0)                   wz = .5 * (1+delta/dz0) * (dz0+delta);
			else if (delta < dz0)                 wz = .5 * delta/dz0 * (2*dz0-delta) + dz0/2.;
			else                                  wz = dz0;

			for (ii=0; ii<=nif[i]; ii++) {

				x0 = (i0f[i]+ii) * dx0;
				i0 = (i0f[i]+ii+Nx0) % Nx0;

				if      ((delta = x-dx1/2.-x0) > dx0) wx = 0; // not likely
				else if (delta > 0)                   wx = .5 * (1-delta/dx0) * (dx0-delta);
				else if (delta > -dx0)                wx = .5 * (-delta/dx0) * (2*dx0+delta) + dx0/2.;
				else if ((delta = x+dx1/2.-x0) <-dx0) wx = 0; // not likely
				else if (delta < 0)                   wx = .5 * (1+delta/dx0) * (dx0+delta);
				else if (delta < dx0)                 wx = .5 * delta/dx0 * (2*dx0-delta) + dx0/2.;
				else                                  wx = dx0;

				weight += wx * wz;
				integr += src.id(i0,j0,k0) * wx * wz;
		}}

		dst.id(i,j1,k) = integr / weight;
	}}
}


void Interp::layerTriFlt(int j0, int j1)
/* planar filter in X & Z direction, defined as weighed average of 9 neighbour nodes */
/* from j0 layer of src to j1 layer of dst, src and dst must have the same planar grid */
{
	if (Nx0 != Nx1 || Nz0 != Nz1)
	{ cout << "Invalid mesh for TriFlt !" << endl; exit(0); }

	int i, j, k, idx, im0, ip0, km0, kp0;
	int im, ip, km, kp, imkm, ipkp, imkp, ipkm;
	double *q0 = src.lyrGet(j0), *q1 = dst.lyrGet(j1);

	for (k=0; k<Nz0; k++) { km0 = ms0.kma[k]; kp0 = ms0.kpa[k];
	for (i=0; i<Nx0; i++) { im0 = ms0.ima[i]; ip0 = ms0.ipa[i];

		idx = Nx0 * k + i;
		im = Nx0 * k + im0;
		ip = Nx0 * k + ip0;
		km = Nx0 * km0 + i;
		kp = Nx0 * kp0 + i;
		imkm = Nx0 * km0 + im0;
		ipkp = Nx0 * kp0 + ip0;
		imkp = Nx0 * kp0 + im0;
		ipkm = Nx0 * km0 + ip0;

		// q1[idx] = 0.25 * q0[idx] + \
		// 	0.125 * (q0[im] + q0[ip] + q0[km] + q0[kp]) + \
		// 	0.0625 * (q0[imkm] + q0[ipkp] + q0[imkp] + q0[ipkm]);
		q1[idx] = 4./9. * q0[idx] + \
			1./9. * (q0[im] + q0[ip] + q0[km] + q0[kp]) + \
			1./36. * (q0[imkm] + q0[ipkp] + q0[imkp] + q0[ipkm]);
	}}
}


void Interp::layerY(int j1, char stgtyp)
/* linear interpolation in Y direction */
/* from src to dst, whole bulk interped */
/* out-of-range points will be linearly extrapolated */
{
	if (Nx0 != Nx1 || Nz0 != Nz1)
	{ cout << "Invalid mesh for y interpolation !" << endl; exit(0); }

	int i, k, j0;
	double y, y1, y2;

	if (stgtyp == 'U')
	{	j0 = j0u[j1];
		y = ms1.yc[j1];
		y1 = ms0.yc[j0];
		y2 = ms0.yc[j0+1];}
	else if (stgtyp == 'V')
	{	j0 = j0v[j1];
		y = ms1.y[j1];
		y1 = ms0.y[j0];
		y2 = ms0.y[j0+1];	}
	else if (stgtyp == 'X')
	{	j0 = j0x[j1];
		y = ms1.yc[j1];
		y1 = ms0.y[j0];
		y2 = ms0.y[j0+1];	}

	for (k=0; k<Nz1; k++) {
	for (i=0; i<Nx1; i++) {
		dst.id(i,j1,k) = 1./(y2-y1) * (
			(y2-y) * src.id(i,j0,k)
		+	(y-y1) * src.id(i,j0+1,k)	);
	}}
}


void Interp::bulkInterp(char stgtyp)
/* combination of planar and wall-normal interpolation, planar first */
{
	int j;
	Mesh ms(Nx1,Ny0,Nz1,Lx1,Ly0,Lz1);
	for (j=0; j<=Ny0; j++) {
		ms.y [j] = ms0.y [j];
		ms.yc[j] = ms0.yc[j];
	}

	Scla mid(ms);
	Interp int1(src,mid), int2(mid,dst);

	for (j=0; j<=Ny0; j++) int1.layerPrdLin(j, j);
	for (j=0; j<=Ny1; j++) int2.layerY(j, stgtyp);

	ms.freeall();
}

void Interp::bulkFilter(char stgtyp)
/* combination of wall-normal interpolation and planar filter, wall-normal first */
{
	int j;
	Mesh ms(Nx0,Ny1,Nz0,Lx0,Ly1,Lz0);
	for (j=0; j<=Ny1; j++) {
		ms.y [j] = ms1.y [j];
		ms.yc[j] = ms1.yc[j];
	}

	Scla mid(ms);
	Interp int1(src,mid), int2(mid,dst);

	for (j=0; j<=Ny1; j++) {
		int1.layerY(j, stgtyp);
		int2.layerPrdLin(j, j); // should change to filter some day
	}

	ms.freeall();
}


// Scla& Scla::interpolate(Scla &src)
// /* initiate flow field from given fields, interpolated to the current grid */
// {
// 	int i, j, k, idx, j0, *i0 = new int [Nx], *k0 = new int [Nz];
// 	Mesh &ms0 = src.meshGet();
// 	int Nx0 = ms0.Nx, Ny0 = ms0.Ny, Nz0 = ms0.Nz, Nxz0 = ms0.Nxz;
// 	double *yc0 = ms0.yc, rescale = (ms0.Lx * ms0.Lz) / (Lx * Lz); // rescale = (alfa * beta) / (alfa0 * beta0)

// 	// nearest interpolation for x and z directions in Fourier space
// 	// find out nearest point for every wavenumber k_z
// 	for (k=0; k<Nz; k++) {	k0[k] = 0;
// 	for (int k1=1; k1<Nz0; k1++) {
// 		if (  fabs(kz(k) - ms0.kz(k1)   )
// 			< fabs(kz(k) - ms0.kz(k0[k]))	)	k0[k] = k1;

// 	}}
// 	// find out nearesr point for every wavenumber k_x
// 	for (i=0; i<Nx; i++) {	i0[i] = 0;
// 	for (int i1=1; i1<Nx0; i1++) {
// 		if (  fabs(kx(i) - ms0.kx(i1)   )
// 			< fabs(kx(i) - ms0.kx(i0[i]))	)	i0[i] = i1;
// 	}}

// 	src.fft();

// 	// linear interpolations for y direction
// 	for (j=0; j<=Ny; j++) {
// 		// find out the range (j0-1 ~ j0) where j fall
// 		for (j0=0; j0<=Ny0; j0++) { if (yc0[j0] > yc[j]) break; }
// 		// if (n == 1) { j0 = j==0 ? 0 : (j==1 ? 2 : j0); }

// 		for (k=0; k<Nz; k++) {
// 		for (i=0; i<(int)(Nx/2+1); i++) {
// 			// values at both ends are extended for out-of-range points
// 			if (j0 == 0) {
// 				this->idf(2*i,  j,k) = src.idf(2*i0[i],  j0,k0[k]);	// real part
// 				this->idf(2*i+1,j,k) = src.idf(2*i0[i]+1,j0,k0[k]);	// imaginary part
// 			}
// 			else if (j0 == Ny0+1) {
// 				this->idf(2*i,  j,k) = src.idf(2*i0[i],  j0-1,k0[k]);
// 				this->idf(2*i+1,j,k) = src.idf(2*i0[i]+1,j0-1,k0[k]);
// 			}
// 			// linear interpolation for in-range points
// 			else {
// 				this->idf(2*i,j,k) =
// 					(yc0[j0] - yc[j]) / ( yc0[j0] - yc0[j0-1] ) * src.idf(2*i0[i],j0-1,k0[k])
// 				+	(yc[j]-yc0[j0-1]) / ( yc0[j0] - yc0[j0-1] ) * src.idf(2*i0[i],j0,  k0[k]);

// 				this->idf(2*i+1,j,k) =
// 					(yc0[j0] - yc[j]) / ( yc0[j0] - yc0[j0-1] ) * src.idf(2*i0[i]+1,j0-1,k0[k])
// 				+	(yc[j]-yc0[j0-1]) / ( yc0[j0] - yc0[j0-1] ) * src.idf(2*i0[i]+1,j0,  k0[k]);
// 			}
// 		}}
// 	}
// 	this->ifft();
// 	this->bulkMlt(1.0 * Nxz / Nxz0);

// 	delete [] i0; delete [] k0;
// 	return *this;
// }



// this function does not belong to the Matrix class
double linearInterp(int n, double *x0, double *y0, double x)
{
	int i;
	for (i=0; i<n; i++) { if (x0[i] > x) break; }
	return	(i==0) ? y0[0]
		:	(i==n) ? y0[n-1]
		:	(	( x0[i] - x ) * y0[i-1]
			+	( x - x0[i-1] ) * y0[i]	) / ( x0[i] - x0[i-1] );
}



// # define DEBUG
# ifdef DEBUG

int main()
{
	Scla src(Mesh(4, 0, 4, 8., 0., 4.));
	Scla dst(Mesh(4, 0, 4, 8., 0., 4.));
	Interp interp(src, dst);

	src.id(0,0,0) = 1; src.id(1,0,0) = 2; src.id(2,0,0) = 3; src.id(3,0,0) = 4;
	src.id(0,0,1) = 5; src.id(1,0,1) = 6; src.id(2,0,1) = 7; src.id(3,0,1) = 8;
	src.id(0,0,2) = 9; src.id(1,0,2) =10; src.id(2,0,2) =11; src.id(3,0,2) =12;
	src.id(0,0,3) =13; src.id(1,0,3) =14; src.id(2,0,3) =15; src.id(3,0,3) =16;

	src.debug_AsciiOutput("", "src", 0, 1);

	interp.layerTriFlt(0, 0);

	dst.debug_AsciiOutput("", "dst", 0, 1);
}

# endif


// make file example
// ##### path names #####
// # header files, source files, middle files and data files
// # relative to Makefile

// IDIR = include
// SDIR = src
// ODIR = obj
// DATADIR = ../data


// ##### compile & link options #####

// OUTPUT = test

// CC = icc
// LIBS = -lfftw3 -lm -qopenmp
// CFLAGS = -I $(IDIR) $(LIBS)


// ##### file lists #####

// _DEPS = Basic.h Interp.h
// DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

// _OBJS = Bulk.o Mesh.o Scla.o Vctr.o Interp.o
// OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))


// #### targets & commands #####
// # note: $@ target, $^ all dependents, $< the first dependent

// $(OUTPUT): $(OBJS)
// 	$(CC) -o $@ $^ $(CFLAGS)

// $(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
// 	$(CC) -c -o $@ $< $(CFLAGS)

// dir:
// 	mkdir -p $(DATADIR)/fielddata
// 	mkdir -p $(DATADIR)/probedata
// 	mkdir -p $(DATADIR)/statdata
// 	mkdir -p $(DATADIR)/postdata

// dirclean:
// 	rm -r $(DATADIR)/fielddata
// 	rm -r $(DATADIR)/probedata
// 	rm -r $(DATADIR)/statdata
// 	rm -r $(DATADIR)/postdata
// 	make dir

// .PHONY: remake
// remake:
// # 	make clean
// 	make
// 	mv $(OUTPUT) $(DATADIR)/

// .PHONY: clean
// clean:
// 	rm -f $(OUTPUT) *~ $(ODIR)/*.o $(INCDIR)/*~

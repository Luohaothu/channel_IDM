# include <iostream>
# include <stdlib.h>
# include <cmath>

# include "Statis.h"

using namespace std;


Statis::Statis(const Mesh &mesh): Mesh(mesh)
{	
	Um = new double [Ny+1];
	Vm = new double [Ny+1];
	Wm = new double [Ny+1];
	Pm = new double [Ny+1];
	R11= new double [Ny+1];
	R22= new double [Ny+1];
	R33= new double [Ny+1];
	R12= new double [Ny+1];
	R23= new double [Ny+1];
	R13= new double [Ny+1];
	Rpu= new double [Ny+1];
	Rpv= new double [Ny+1];
	Rpw= new double [Ny+1];
	Rpp= new double [Ny+1];
	Num= new double [Ny+1];
}

Statis::~Statis()
{
	delete [] Um; delete [] Vm; delete [] Wm; delete [] Pm;
	delete [] R11;delete [] R22;delete [] R33;
	delete [] R12;delete [] R23;delete [] R13;
	delete [] Rpu;delete [] Rpv;delete [] Rpw;delete [] Rpp;
	delete [] Num;
}

double Statis::checkDiv(const Vctr &U)
{
	int i, j, k; double divmax = -1.;
	const Mesh &ms = U.meshGet();

	for (j=1; j<ms.Ny; j++) {
	for (k=0; k<ms.Nz; k++) {
	for (i=0; i<ms.Nx; i++) {
		div = fabs(U.divergence(i, j, k));
		if ( div > divmax || divmax < 0 ) {
			divmax = div;
			divpos[0] = i;
			divpos[1] = j;
			divpos[2] = k;
		}
	}}}
	return ( div = divmax );
}

double Statis::checkCFL(const Vctr &U, double dt)
{
	int i, j, k; double cflmax = -1.;
	const Mesh &ms = U.meshGet();

	for (j=1; j<ms.Ny; j++) {
	for (k=0; k<ms.Nz; k++) {
	for (i=0; i<ms.Nx; i++) {
		cfl = U.module(i, j, k) * dt;
		if ( cfl > cflmax || cflmax < 0 ) {
			cflmax = cfl;
			cflpos[0] = i;
			cflpos[1] = j;
			cflpos[2] = k;
		}
	}}}
	return ( cfl = cflmax );
}

double Statis::checkTaub(const Vctr &U, double Re)
/* calculate 3 components of total stress acting on the boundary: tau_2i = nu d<u_i>/dy - <v'u_i'> */
/* tau_2i = 2 nu S_2i - <v'u_i'>; 2*S_2i = du_i/dy + dv/dx_i = d<u_i>/dy (applying mass conservation) */
{
	Mesh ms(Nx,0,Nz,Lx,0,Lz);
	Scla ul1(ms), ul2(ms), ul3(ms), ul4(ms),
	     vl1(ms), vl2(ms), vl3(ms), vl4(ms),
	     wl1(ms), wl2(ms), wl3(ms), wl4(ms);
	double sm21, sm22, sm23, rm21, rm22, rm23;

	U.layer2CC(ul1[0], vl1[0], wl1[0], 0);
	U.layer2CC(ul2[0], vl2[0], wl2[0], 1);
	U.layer2CC(ul3[0], vl3[0], wl3[0], Ny-1);
	U.layer2CC(ul4[0], vl4[0], wl4[0], Ny);

	sm21 = (ul2.av() - ul1.av()) / dy[1] - (ul4.av() - ul3.av()) / dy[Ny-1]; // d<u>/dy (0.5 compensated by dy)
	sm22 = (vl2.av() - vl1.av()) / dy[1] + (vl4.av() - vl3.av()) / dy[Ny-1]; // d<v>/dy
	sm23 = (wl2.av() - wl1.av()) / dy[1] - (wl4.av() - wl3.av()) / dy[Ny-1]; // d<w>/dy

	(vl2 = vl1) += -vl1.av(); // v' bottom
	(vl3 = vl4) += -vl4.av(); // v' top
	rm21 = ( (ul1 += -ul1.av()) *= vl2 ).av() - ( (ul4 += -ul4.av()) *= vl3 ).av(); // 2<u'v'>
	rm22 = ( (vl1 += -vl1.av()) *= vl2 ).av() + ( (vl4 += -vl4.av()) *= vl3 ).av(); // 2<v'v'>
	rm23 = ( (wl1 += -wl1.av()) *= vl2 ).av() - ( (wl4 += -wl4.av()) *= vl3 ).av(); // 2<w'v'>

	taub[0] = sm21 / Re - rm21 / 2.;
	taub[1] = sm22 / Re - rm22 / 2.;
	taub[2] = sm23 / Re - rm23 / 2.;

	ms.freeall();
	return taub[0];

	// taub[0] = taub[1] = taub[2] = 0;
	// double u1, u2, u3, u4, v1, v2, v3, v4, w1, w2, w3, w4;
	// double um1 = ul1.layerMean(), um4 = ul4.layerMean();
	// double vm1 = vl1.layerMean(), vm4 = vl4.layerMean();
	// double wm1 = wl1.layerMean(), wm4 = wl4.layerMean();
	// for (k=0; k<Nz; k++) {
	// for (i=0; i<Nx; i++) {
	// 	u1 = ul1.id(i,0,k); u2 = ul2.id(i,0,k); u3 = ul3.id(i,0,k); u4 = ul4.id(i,0,k);
	// 	v1 = vl1.id(i,0,k); v2 = vl2.id(i,0,k); v3 = vl3.id(i,0,k); v4 = vl4.id(i,0,k);
	// 	w1 = wl1.id(i,0,k); w2 = wl2.id(i,0,k); w3 = wl3.id(i,0,k); w4 = wl4.id(i,0,k);

	// 	taub[0] += ( 1./Re * ( (u2-u1)/dy[1] - (u4-u3)/dy[Ny-1] ) - .5 * ( (u1-um1)*(v1-vm1) - (u4-um4)*(v4-vm4) ) ) / Nxz;
	// 	taub[1] += ( 1./Re * ( (v2-v1)/dy[1] + (v4-v3)/dy[Ny-1] ) - .5 * ( (v1-vm1)*(v1-vm1) + (v4-vm4)*(v4-vm4) ) ) / Nxz;
	// 	taub[2] += ( 1./Re * ( (w2-w1)/dy[1] - (w4-w3)/dy[Ny-1] ) - .5 * ( (w1-wm1)*(v1-vm1) - (w4-wm4)*(v4-vm4) ) ) / Nxz;
	// }}
}

double Statis::checkEner(const Vctr &U, const Scla &P, const Scla &NU)
/* calculate Reynolds stresses, pressure-velocity correlations, and the total fluctuation energy */
{
	Mesh ms(Nx,0,Nz,Lx,0,Lz);
	Scla ql(ms), ul(ms), vl(ms), wl(ms), pl(ms);

	velm[0] = 0;
	velm[1] = 0;
	velm[2] = 0;
	ener = 0;

	for (int j=0; j<=Ny; j++) {
		// interpolate to cell centers or real boundaries
		U.layer2CC(ul[0], vl[0], wl[0], j);
		pl = P[j==0 ? 1 : j==Ny ? Ny-1 : j];

		// calculate means and fluctuations at cell centers
		ul += -(Um[j] = ul.av());
		vl += -(Vm[j] = vl.av());
		wl += -(Wm[j] = wl.av());
		pl += -(Pm[j] = pl.av());
		       Num[j] = NU.av(j);

		// calculate Reynolds stress and cross correlations at cell centers
		R11[j] = ( (ql = ul) *= ul ).av();
		R22[j] = ( (ql = vl) *= vl ).av();
		R33[j] = ( (ql = wl) *= wl ).av();
		R12[j] = ( (ql = ul) *= vl ).av();
		R23[j] = ( (ql = vl) *= wl ).av();
		R13[j] = ( (ql = wl) *= ul ).av();
		Rpu[j] = ( (ql = pl) *= ul ).av();
		Rpv[j] = ( (ql = pl) *= vl ).av();
		Rpw[j] = ( (ql = pl) *= wl ).av();
		Rpp[j] = ( (ql = pl) *= pl ).av();

		// wall-normal integration for bulk velocities and turbulent energy
		if (0 < j && j < Ny) {
			double weight = dy[j] / Ly;
			velm[0] += Um[j] * weight;
			velm[1] += Vm[j] * weight;
			velm[2] += Wm[j] * weight;
			ener += .5 * (R11[j] + R22[j] + R33[j]) * weight;
		}
	}

	ms.freeall();
	return ener;
}


void Statis::writeProfile(const char *path, int tstep) const
{
	FILE *fp;
	char str[1024];

	if (tstep < 0)	sprintf(str, "%sPROF.dat", path?path:"");
	else			sprintf(str, "%sPROF%08i.dat", path?path:"", tstep);

	fp = fopen(str, "w");

		fputs("Title = \"Instantaneous profile\"\n", fp);
		fputs("variables = \"y\", \"U\", \"V\", \"W\", \"P\", \"R11\", \"R22\", \"R33\", \"R12\", \"R23\", \"R13\", \"Rpu\", \"Rpv\", \"Rpw\", \"Rpp\", \"NU\"\n", fp);
		sprintf(str, "zone t = \"%i\", i = %i\n", tstep, Ny+1);
		fputs(str, fp);

		for (int j=0; j<=Ny; j++) {
			sprintf(str, "%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
				j==0 ? y[1] : j==Ny ? y[Ny] : yc[j],
				Um[j], Vm[j], Wm[j], Pm[j],
				R11[j], R22[j], R33[j], R12[j], R23[j], R13[j],
				Rpu[j], Rpv[j], Rpw[j], Rpp[j],
				Num[j]	);
			fputs(str, fp);
		}

	fclose(fp);
}


void Statis::writeLogfile(const char *path, int tstep, double time, const double mpg[3]) const
{
	FILE *fp;
	char str[1024];
	long int pos = this->getLogpos(path, tstep);

	sprintf(str, "%sLOG.dat", path?path:"");

	if (pos > 0) {
		fp = fopen(str, "r+");
		fseek(fp, 0, SEEK_END);
		if (pos < ftell(fp)) {
			char *buf = new char [pos / sizeof(char)];
			fseek(fp, 0, SEEK_SET);
			fread(buf, pos, 1, fp);
			fclose(fp);
			fp = fopen(str, "w");
			fwrite(buf, pos, 1, fp);
			delete [] buf;
		}
	}
	else {
		fp = fopen(str, "w");
		fputs("Title = \"Running log\"\n", fp);
		fputs("variables = \"n\", \"t\", \"ener\", \"taub21\", \"taub22\", \"taub23\", \"mpg1\", \"mpg2\", \"mpg3\", \"Um\", \"Vm\", \"Wm\", \"div\", \"divposx\", \"divposy\", \"divposz\", \"cfl\", \"cflposx\", \"cflposy\", \"cflposz\"\n", fp);
		fputs("zone t = \"statis\"\n", fp);
	}

	sprintf(str, "%i\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%i\t%i\t%i\t%.18e\t%i\t%i\t%i\n",
		tstep, time, ener,
		taub[0], taub[1], taub[2],
		mpg [0], mpg [1], mpg [2],
		velm[0], velm[1], velm[2],
		div, divpos[0], divpos[1], divpos[2],
		cfl, cflpos[0], cflpos[1], cflpos[2]	);
	fputs(str, fp);

	fclose(fp);
}


long int Statis::getLogpos(const char *path, int tstep)
/* find the first line whose time step >= tstep, return the beginning position and record the time of this line */
{
	FILE *fp;
	char str[1024];
	long int pos = 0;
	int n = 0;
	double t = 0.0;

	sprintf(str, "%sLOG.dat", path?path:"");

	if ( ( fp = fopen(str, "r") ) ) {
		fgets(str, 1024, fp);	// skip header
		fgets(str, 1024, fp);	// skip header
		fgets(str, 1024, fp);	// skip header

		while ( n < tstep ) {
			pos = ftell(fp);
			if ( ! fgets(str, 1024, fp) ) break; // reach the end of file
			sscanf(str, "%i\t%lf", &n, &t);	// if str is empty line or spaces, n and t will not be assigned
		}
		fclose(fp);
	}
	
	return pos;
}

double Statis::getLogtime(const char *path, int tstep)
{
	FILE *fp;
	char str[1024];
	long int pos = getLogpos(path, tstep);
	int n = 0;
	double t = 0.0;

	sprintf(str, "%sLOG.dat", path?path:"");

	if (pos > 0) {
		fp = fopen(str, "r");
		fseek(fp, pos, SEEK_SET);
		if(fgets(str, 1024, fp)) sscanf(str, "%i\t%lf", &n, &t);
		fclose(fp);
	}

	return t;
}









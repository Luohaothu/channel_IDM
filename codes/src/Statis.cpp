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

double Statis::checkDiv(Vctr &U)
{
	int i, j, k;
	double divmax = -1.0;

	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
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

double Statis::checkCFL(Vctr &U, double dt)
{
	int i, j, k;
	double cflmax = -1.0;

	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		cfl = U.convection(i, j, k) * dt;
		if ( cfl > cflmax || cflmax < 0 ) {
			cflmax = cfl;
			cflpos[0] = i;
			cflpos[1] = j;
			cflpos[2] = k;
		}
	}}}
	return ( cfl = cflmax );
}

double Statis::checkMean(Vctr &U, Scla &P, Scla &NU)
/* calculate mean values at cell centers */
{
	double vm1, vm2;
	Scla &U1 = U.com1, &U2 = U.com2, &U3 = U.com3;

	for (int j=0; j<=Ny; j++) {
		vm1 = j==0 ? U2.layerMean(1) : vm2;
		vm2 = j==Ny ? vm1 : U2.layerMean(j+1);
		Vm[j] = 0.5 * (vm1 + vm2);
		Um[j] = U1.layerMean(j);	// layer mean of U at U grid points is equivalent to that at cell centers
		Wm[j] = U3.layerMean(j);	// layer mean of W at W grid points is equivalent to that at cell centers
		Pm[j] =  P.layerMean(j);	// P stored in cell centers
		Num[j]= NU.layerMean(j);	// NU stored in cell centers
	}
	velm[0] = U1.bulkMeanU();
	velm[1] = U2.bulkMeanV();
	velm[2] = U3.bulkMeanU();
	return velm[0];
}

double Statis::checkTaub(Vctr &U, double Re)
/* calculate three components of total stress acting on the boundary */
/* tau_2i = 2 nu S_2i - <v'u_i'>; S_2i = 0.5 ( du_i/dy + dv/dx_i ) */
/* applying mass conservation, on the boundary tau_2i = nu d<u_i>/dy - <v'u_i'> */
{
	double qm1, qm2;
	Scla q(Mesh(Nx,0,Nz,Lx,0,Lz)), v(q.meshGet());


	U.com2.layerVG2CC(v, 0, 0).layerAdd(-v.layerMean());

	qm2 = U.com1.layerUG2CC(q, 0, 1).layerMean();
	U.com1.layerUG2CC(q, 0, 0).layerAdd(-(qm1 = q.layerMean()));
	taub[0] = 1./Re * (qm2 - qm1) / h[1] - q.layersMlt(v).layerMean();

	qm2 = U.com2.layerVG2CC(q, 0, 1).layerMean();
	U.com2.layerVG2CC(q, 0, 0).layerAdd(-(qm1 = q.layerMean()));
	taub[1] = 1./Re * (qm2 - qm1) / h[1] - q.layersMlt(v).layerMean();

	qm2 = U.com3.layerWG2CC(q, 0, 1).layerMean();
	U.com3.layerWG2CC(q, 0, 0).layerAdd(-(qm1 = q.layerMean()));
	taub[2] = 1./Re * (qm2 - qm1) / h[1] - q.layersMlt(v).layerMean();


	U.com2.layerVG2CC(v, 0, Ny).layerAdd(-v.layerMean());

	qm1 = U.com1.layerUG2CC(q, 0, Ny-1).layerMean();
	U.com1.layerUG2CC(q, 0, Ny).layerAdd(-(qm2 = q.layerMean()));
	taub[0] -= 1./Re * (qm2 - qm1) / h[1] - q.layersMlt(v).layerMean();

	qm1 = U.com2.layerVG2CC(q, 0, Ny-1).layerMean();
	U.com2.layerVG2CC(q, 0, Ny).layerAdd(-(qm2 = q.layerMean()));
	taub[1] -= 1./Re * (qm2 - qm1) / h[1] - q.layersMlt(v).layerMean();

	qm1 = U.com3.layerWG2CC(q, 0, Ny-1).layerMean();
	U.com3.layerWG2CC(q, 0, Ny).layerAdd(-(qm2 = q.layerMean()));
	taub[2] -= 1./Re * (qm2 - qm1) / h[1] - q.layersMlt(v).layerMean();


	taub[0] *= .5;
	taub[1] *= .5;
	taub[2] *= .5;

	q.meshGet().freeall();

	// double taud, tauu, *vel[3] = {Um, Vm, Wm};
	// for (int i=0; i<3; i++) {
	// 	taud = ( vel[i][1]	- vel[i][0]		) / h[1];
	// 	tauu = ( vel[i][Ny] - vel[i][Ny-1]	) / h[Ny];
	// 	taub[i] = 0.5 * ( taud - tauu ) * ( i==1 ? (2.0/Re) : (1.0/Re) );
	// }
	return taub[0];
}

// note: must call checkMean() before this function
double Statis::checkEner(Vctr &U, Scla &P)
/* calculate Reynolds stresses, pressure-velocity correlations, and the total fluctuation energy */
{
	Scla ql(Mesh(Nx,0,Nz,Lx,0,Lz)), \
	ul(ql.meshGet()), vl(ql.meshGet()), wl(ql.meshGet()), pl(ql.meshGet());

	ener = 0.0;

	for (int j=0; j<=Ny; j++) {
		// calculate fluctuations at cell centers
		U.com1.layerUG2CC(ul,0,j).layerAdd(-Um[j]);
		U.com2.layerVG2CC(vl,0,j).layerAdd(-Vm[j]);
		U.com3.layerWG2CC(wl,0,j).layerAdd(-Wm[j]);
		pl.layerCpy(P,0,j).layerAdd(-Pm[j]);

		// calculate Reynolds stress and cross correlations at cell centers
		R11[j] = ql.layerCpy(ul).layersMlt(ul).layerMean();
		R22[j] = ql.layerCpy(vl).layersMlt(vl).layerMean();
		R33[j] = ql.layerCpy(wl).layersMlt(wl).layerMean();
		R12[j] = ql.layerCpy(ul).layersMlt(vl).layerMean();
		R23[j] = ql.layerCpy(vl).layersMlt(wl).layerMean();
		R13[j] = ql.layerCpy(wl).layersMlt(ul).layerMean();

		Rpu[j] = ql.layerCpy(pl).layersMlt(ul).layerMean();
		Rpv[j] = ql.layerCpy(pl).layersMlt(vl).layerMean();
		Rpw[j] = ql.layerCpy(pl).layersMlt(wl).layerMean();
		Rpp[j] = ql.layerCpy(pl).layersMlt(pl).layerMean();

		ener += 0.5 * ( R11[j] + R22[j] + R33[j] ) * ( dy[j] / Ly );
	}

	ql.meshGet().freeall();
	return ener;
}


void Statis::writeProfile(char *path, int tstep)
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
				yc[j],
				Um[j], Vm[j], Wm[j], Pm[j],
				R11[j], R22[j], R33[j], R12[j], R23[j], R13[j],
				Rpu[j], Rpv[j], Rpw[j], Rpp[j],
				Num[j]	);
			fputs(str, fp);
		}

	fclose(fp);
}





long int Statis::getLogpos(char *path, int tstep)
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

double Statis::getLogtime(char *path, int tstep)
{
	FILE *fp;
	char str[1024];
	long int pos = this->getLogpos(path, tstep);
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

void Statis::writeLogfile(char *path, int tstep, double time)
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









# include <iostream>


# include "Statis.h"

using namespace std;


Statis::Statis(int *dim):
Nx(dim[0]), Ny(dim[1]), Nz(dim[2]), Nxz(dim[0]*dim[2])
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
}

Statis::~Statis()
{
	delete [] Um; delete [] Vm; delete [] Wm; delete [] Pm;
	delete [] R11;delete [] R22;delete [] R33;
	delete [] R12;delete [] R23;delete [] R13;
	delete [] Rpu;delete [] Rpv;delete [] Rpw;delete [] Rpp;
}

double Statis::checkDiv(double **U, class Mesh *pmesh)
{
	int i, j, k;
	double divmax = -1.0, *u = U[0], *v = U[1], *w = U[2];

	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		div = abs( pmesh->divergence(u, v, w, i, j, k) );
		if ( div > divmax || divmax < 0 ) {
			divmax = div;
			divpos[0] = i;
			divpos[1] = j;
			divpos[2] = k;

		}
	}}}
	return ( div = divmax );
}

double Statis::checkCFL(double **U, class Mesh *pmesh, double dt)
{
	int i, j, k;
	double cflmax = -1.0, *u = U[0], *v = U[1], *w = U[2];

	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		cfl = pmesh->convection(u, v, w, i, j, k) * dt;
		if ( cfl > cflmax || cflmax < 0 ) {
			cflmax = cfl;
			cflpos[0] = i;
			cflpos[1] = j;
			cflpos[2] = k;
		}
	}}}
	return ( cfl = cflmax );
}

double Statis::checkMean(double **UP, class Mesh *pmesh)
/* calculate mean values at cell centers */
{
	double vm1, vm2, *u = UP[0], *v = UP[1], *w = UP[2], *p = UP[3];

	for (int j=0; j<=Ny; j++) {
		vm1 = j==0 ? pmesh->layerMean(v, 1) : vm2;
		vm2 = j==Ny ? vm1 : pmesh->layerMean(v, j+1);
		Vm[j] = 0.5 * (vm1 + vm2);
		Um[j] = pmesh->layerMean(u, j);	// layer mean of U at U grid points is equivalent to that at cell centers
		Wm[j] = pmesh->layerMean(w, j);	// layer mean of W at W grid points is equivalent to that at cell centers
		Pm[j] = pmesh->layerMean(p, j);
	}
	velm[0] = pmesh->bulkMeanU(u);
	velm[1] = pmesh->bulkMeanV(v);
	velm[2] = pmesh->bulkMeanU(w);
	return velm[0];
}

double Statis::checkTauw(class Mesh *pmesh)
/* calculate three components of wall stress */
{	// note: must call checkMean() before this function
	double taud, tauu, *vel[3] = {Um, Vm, Wm};	for (int i=0; i<3; i++) {
		taud = ( vel[i][1]	- vel[i][0]		) / pmesh->h[1];
		tauu = ( vel[i][Ny] - vel[i][Ny-1]	) / pmesh->h[Ny];
		tauw[i] = 0.5 * ( taud - tauu );
	}
	return tauw[0];
}

double Statis::checkEner(double **UP, class Mesh *pmesh, class Field *pfield)
/* calculate Reynolds stresses, pressure-velocity correlations, and the total fluctuation energy */
{	// note: must call checkMean() before this function
	int i, j, k, idx, ip, kp;
	double *u = UP[0], *v = UP[1], *w = UP[2], *p = UP[3];
	double *ul = new double [Nxz];
	double *vl = new double [Nxz];
	double *wl = new double [Nxz];
	double *pl = new double [Nxz];
	// double vm1, vm2;
	double *vl1 = vl, *vl2 = new double [Nxz];
	double *ql = new double [Nxz];
	
	ener = 0.0;

	for (j=0; j<=Ny; j++) {
		// calculate fluctuations at cell centers
		vl1 = j==0 ? pfield->layerCopy(vl2,v,0,1) : vl2;
		vl2 = j==Ny ? pfield->layerCopy(vl,vl1) : pfield->layerCopy(vl,v,0,j+1);
		vl = pfield->layerMult( pfield->layersAdd(vl1,vl2), 0.5 );
		pfield->layerAdd( vl, -Vm[j] );
		pfield->layerAdd( pmesh->layerCenterU(ul,u,0,j), -Um[j] );
		pfield->layerAdd( pmesh->layerCenterW(wl,w,0,j), -Wm[j] );
		pfield->layerAdd( pfield->layerCopy(pl,p,0,j), -Pm[j] );

		// calculate Reynolds stress and cross correlations at cell centers
		R11[j] = pmesh->layerMean( pfield->layersMult( pfield->layerCopy(ql,ul) ) );
		R22[j] = pmesh->layerMean( pfield->layersMult( pfield->layerCopy(ql,vl) ) );
		R33[j] = pmesh->layerMean( pfield->layersMult( pfield->layerCopy(ql,wl) ) );
		R12[j] = pmesh->layerMean( pfield->layersMult( pfield->layerCopy(ql,ul), vl ) );
		R23[j] = pmesh->layerMean( pfield->layersMult( pfield->layerCopy(ql,vl), wl ) );
		R13[j] = pmesh->layerMean( pfield->layersMult( pfield->layerCopy(ql,wl), ul ) );

		Rpu[j] = pmesh->layerMean( pfield->layersMult( pfield->layerCopy(ql,pl), ul ) );
		Rpv[j] = pmesh->layerMean( pfield->layersMult( pfield->layerCopy(ql,pl), vl ) );
		Rpw[j] = pmesh->layerMean( pfield->layersMult( pfield->layerCopy(ql,pl), wl ) );
		Rpp[j] = pmesh->layerMean( pfield->layersMult( pfield->layerCopy(ql,pl) ) );

		ener += 0.5 * ( R11[j] + R22[j] + R33[j] ) * ( pmesh->dy[j] / pmesh->Ly );
	}
	return ener;
}


void Statis::writeProfile(char *path, int tstep, double *yc)
{
	FILE *fp;
	char str[1024];

	if (tstep < 0)	sprintf(str, "%sPROF.dat", path?path:"");
	else			sprintf(str, "%sPROF%08i.dat", path?path:"", tstep);

	fp = fopen(str, "w");

		fputs("Title = \"Instantaneous profile\"\n", fp);
		fputs("variables = \"y\", \"U\", \"V\", \"W\", \"P\", \"R11\", \"R22\", \"R33\", \"R12\", \"R23\", \"R13\", \"Rpu\", \"Rpv\", \"Rpw\", \"Rpp\"\n", fp);
		sprintf(str, "zone t = \"%i\", i = %i\n", tstep, Ny+1);
		fputs(str, fp);

		for (int j=0; j<=Ny; j++) {
			sprintf(str, "%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
				yc[j],
				Um[j], Vm[j], Wm[j], Pm[j],
				R11[j], R22[j], R33[j], R12[j], R23[j], R13[j],
				Rpu[j], Rpv[j], Rpw[j], Rpp[j]	);
			fputs(str, fp);
		}

	fclose(fp);
}


void Statis::writeLogfile(char *path, int tstep, double time)
{
	FILE *fp;
	char str[1024];
	bool fileexist = true;

	sprintf(str, "%sLOG.dat", path?path:"");

	if ( ( fp = fopen(str, "r") ) ) fclose(fp);
	else fileexist = false;

	fp = fopen(str, "a");
	
		if (! fileexist) {
			fputs("Title = \"Running log\"\n", fp);
			fputs("variables = \"n\", \"t\", \"ener\", \"tauw21\", \"tauw22\", \"tauw23\", \"Um\", \"Vm\", \"Wm\", \"div\", \"divposx\", \"divposy\", \"divposz\", \"cfl\", \"cflposx\", \"cflposy\", \"cflposz\"\n", fp);
			fputs("zone t = \"statis\"\n", fp);
		}

		sprintf(str, "%i\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%i\t%i\t%i\t%.18e\t%i\t%i\t%i\n",
			tstep, time, ener,
			tauw[0], tauw[1], tauw[2],
			velm[0], velm[1], velm[2],
			div, divpos[0], divpos[1], divpos[2],
			cfl, cflpos[0], cflpos[1], cflpos[2]	);
		fputs(str, fp);

	fclose(fp);
}










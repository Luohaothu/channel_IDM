# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <cmath>


# include "Basic.h"

using namespace std;



Mesh::Mesh(int n1, int n2, int n3, double l1, double l2, double l3, double dy_min, char* path):
/* allocate memory for pointers */
Nx(n1), Ny(n2), Nz(n3), Nxz(n1 * n3),
Lx(l1), Ly(l2), Lz(l3), vol(Lx * Ly * Lz),
dx((double) Lx / Nx), dx2(dx * dx),
dz((double) Lz / Nz), dz2(dz * dz)
{
	y = new double [Ny+1];	// Ny layers (including 2 walls) + 1 redundant layer (index 0)
	yc= new double [Ny+1];	// Ny+1 layers (including 2 walls)

	dy= new double [Ny+1];	// Ny-1 boxes + 2 redundant (index 0 and Ny)
	h = new double [Ny+1];	// Ny layers + 1 redundant layer (index 0)

	dvol = new double [Ny+1];

	hm = new double [Ny];
	hc = new double [Ny];
	hp = new double [Ny];

	dym = new double [Ny];
	dyc = new double [Ny];
	dyp = new double [Ny];

	pmj = new double [Ny];
	pcj = new double [Ny];
	ppj = new double [Ny];

	ak1 = new double [Nx];
	ak3 = new double [Nz];

	kpa = new int [Nz];
	kma = new int [Nz];

	ipa = new int [Nx];
	ima = new int [Nx];

	int i, k;
	// initiate wavenumbers
	for (i=0; i<Nx; i++) ak1[i] = 2.0/dx2 * ( 1.0 - cos(kx(i) * dx) );
	for (k=0; k<Nz; k++) ak3[k] = 2.0/dz2 * ( 1.0 - cos(kz(k) * dz) );
	// initiate indexes
	for (k=0; k<Nz; k++) {
		kpa[k] = (k+1) % Nz;
		kma[k] = (k-1+Nz) % Nz;
	}
	for (i=0; i<Nx; i++) {
		ipa[i] = (i+1) % Nx;
		ima[i] = (i-1+Nx) % Nx;
	}
	// initiate y mesh
	if (dy_min > 0){ this->getYmesh(dy_min);this->initYmesh(); }
	else if (path) { this->getYmesh(path);  this->initYmesh(); }
}

void Mesh::freeall()
{
	delete [] y; delete [] yc;
	delete [] dy; delete [] h; delete [] dvol;
	delete [] hm; delete [] hc; delete [] hp;
	delete [] dym; delete [] dyc; delete [] dyp;
	delete [] pmj; delete [] pcj; delete [] ppj;
	delete [] ak1; delete [] ak3;
	delete [] kpa; delete [] kma;
	delete [] ipa; delete [] ima;
}


void Mesh::initYmesh()
{
	int i, j, k;

	// check
	if ( fabs((y[Ny]-y[1]) - Ly) > 1e-10 )
		{ cout << endl << "No valid Y mesh provided !" << endl; exit(0); }

	// y coordinates
	y[0] = 0;
	yc[0] = y[1];
	yc[Ny] = y[Ny];
	for (j=1; j<Ny; j++) yc[j] = 0.5 * ( y[j] + y[j+1] );

	// y intervals
	for (j=1; j<Ny; j++) { dy[j] = y[j+1] - y[j]; dvol[j] = dy[j] * dx * dz; }
	for (j=2; j<Ny; j++) { h[j] = (y[j+1] - y[j-1]) / 2.0; }
	dy[0] = dvol[0] = 0.0;
	dy[Ny]= dvol[Ny]= 0.0;
	h[0] = 0.0;
	h[1] = dy[1] / 2.0;
	h[Ny] = dy[Ny-1] / 2.0;

	// coefficients for wall-normal 2nd-order derivative
	for (j=1; j<Ny; j++) {	// for U & W in laplacian operator, boundaries not excluded
		hp[j] = 1.0 / dy[j] / h[j+1];
		hc[j] = 1.0 / dy[j] * ( 1.0/h[j+1] + 1.0/h[j] );	// need to multiply by -1 in use
		hm[j] = 1.0 / dy[j] / h[j];
	}
	for (j=2; j<Ny; j++) {	// for V in laplacian operator, boundaries not excluded
		dyp[j] = 1.0 / h[j] / dy[j];
		dyc[j] = 1.0 / h[j] * ( 1.0/dy[j] + 1.0/dy[j-1] );	// need to multiply by -1 in use
		dym[j] = 1.0 / h[j] / dy[j-1];
	}
	for (j=1; j<Ny; j++) {	// for P in divergence-gradient operator, boundaries excluded
		if (j>1 && j<Ny-1) {
			ppj[j] = hp[j];
			pcj[j] = - hc[j];								// NO need to multiply by -1 in use
			pmj[j] = hm[j];
		}
		else if (j==1) {
			ppj[j] = hp[1];
			pcj[j] = - 1.0 / dy[1] / h[2];
			pmj[j] = 0.0;
		}
		else if (j==Ny-1) {
			ppj[j] = 0.0;
			pcj[j] = - 1.0 / dy[Ny-1] / h[Ny-1];
			pmj[j] = hm[Ny-1];
		}
	}
}

void Mesh::writeYmesh(char *path)
{
	char str[1024];
	FILE *fp = fopen(strcat(strcpy(str, path), "CHANNEL.GRD"), "w");
	for (int j=0; j<=Ny; j++) {
		sprintf(str, "%.18e\n", y[j]);
		fputs(str, fp);	}
	fclose(fp);
	cout << "\ndy_min = " << y[2]-y[1] << ", dy_max = " << y[Ny/2+1]-y[Ny/2] << endl;
}



double hyptan(int j, double gamma, int ny, double ly)
/* generate y grids with hyperbolic tangent distribution */
{
	double delta = ly / 2.0;
	double ytild = 2.0 * (j-1) / (ny-1) - 1.0;
	return 1.0 + delta * tanh(gamma * delta * ytild) / tanh(gamma * delta);
}
void Mesh::getYmesh(double dy_min)
/* generate hyperbolic tangent grid (Newton iteration to determine parameter)
target equation: F(gamma) = hyptan(2, gamma, Ny, Ly) - hyptan(1, gamma, Ny, l2) - dy_min = 0 */
{
	double F, gamma = 1.0, dgamma = 0.1, grad;
	double F0 = hyptan(2, gamma-dgamma, Ny, Ly) - hyptan(1, gamma-dgamma, Ny, Ly) - dy_min;
	// determine hyperbolic tangent parameter
	for (int j=0; j<100; j++) {
		F = hyptan(2, gamma, Ny, Ly) - hyptan(1, gamma, Ny, Ly) - dy_min;
		grad = (F-F0) / dgamma;	if (! grad) break;
		dgamma = - F / grad;	if (! dgamma) break;
		gamma += dgamma;		if (gamma <= 0) { cout << "Mesh error: dy_min too large !" << endl; exit(0); }
		F0 = F;
	}
	// generate grid
	y[0] = 0;
	for (int j=1; j<=Ny; j++)
		y[j] = hyptan(j, gamma, Ny, Ly);
}

void Mesh::getYmesh(char *path)
/* read grid from file at specified path */
{
	char str[1024];
	FILE *fp = fopen(strcat(strcpy(str, path), "CHANNEL.GRD"), "r");
	for (int j=0; j<=Ny; j++)
		sscanf(fgets(str, 1024, fp), "%le", & y[j]);
	fclose(fp);
}


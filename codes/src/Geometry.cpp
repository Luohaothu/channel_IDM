#include "Basic.h"

using namespace std;

static double hyptan(int i, double gma, int N, double L);
static double get_dymin(double gma, int N, double L);
static double newton_iter(double tgt, double (*pfun)(double, int, double), int N, double L);



Geometry::Geometry(int Nx, int Ny, int Nz, double Lx, double Ly, double Lz):
Nx(Nx),
Ny(Ny),
Nz(Nz),
Lx(Lx),
Ly(Ly),
Lz(Lz),
Nxc((Nx-1)/2+1),
Nzc((Nz-1)/2+1),
Nxr(Nxc*2),
Nzr(Nzc*2)
{
	x = new double[Nx+1]; xc = new double[Nx+1];
	y = new double[Ny+1]; yc = new double[Ny+1];
	z = new double[Nz+1]; zc = new double[Nz+1];

	dx = new double[Nx+1]; hx = new double[Nx+1];
	dy = new double[Ny+1]; hy = new double[Ny+1];
	dz = new double[Nz+1]; hz = new double[Nz+1];

	kx = new double[Nx]; kx2 = new double[Nx];
	kz = new double[Nz]; kz2 = new double[Nz];

	ima = new int[Nx+1]; ipa = new int[Nx+1];
	jma = new int[Ny+1]; jpa = new int[Ny+1];
	kma = new int[Nz+1]; kpa = new int[Nz+1];
}

Geometry::~Geometry()
{
	delete[] x; delete[] xc;
	delete[] y; delete[] yc;
	delete[] z; delete[] zc;

	delete[] dx; delete[] hx;
	delete[] dy; delete[] hy;
	delete[] dz; delete[] hz;

	delete[] kx; delete[] kx2;
	delete[] kz; delete[] kz2;

	delete[] ima; delete[] ipa;
	delete[] jma; delete[] jpa;
	delete[] kma; delete[] kpa;
}

Geometry::Geometry(const Geometry &geo):
Geometry(geo.Nx,geo.Ny,geo.Nz,geo.Lx,geo.Ly,geo.Lz)
{
	DeepCopy(geo);
}

Geometry& Geometry::operator=(const Geometry &geo)
{
	if (&geo == this) return *this;
	if (Nx != geo.Nx || Ny != geo.Ny || Nz != geo.Nz ||
		Lx != geo.Lx || Ly != geo.Ly || Lz != geo.Lz)
	{
		cout << "Geometry sizes do not match !" << endl;
		exit(0);
	}
	DeepCopy(geo);
	return *this;
}

void Geometry::DeepCopy(const Geometry &geo)
{
	for (int n=0; n<=max(max(Nx,Ny),Nz); n++) {
		if (n<=Nx) { x[n] = geo.x[n]; xc[n] = geo.xc[n]; }
		if (n<=Ny) { y[n] = geo.y[n]; yc[n] = geo.yc[n]; }
		if (n<=Nz) { z[n] = geo.z[n]; zc[n] = geo.zc[n]; }

		if (n<=Nx) { dx[n] = geo.dx[n]; hx[n] = geo.hx[n]; }
		if (n<=Ny) { dy[n] = geo.dy[n]; hy[n] = geo.hy[n]; }
		if (n<=Nz) { dz[n] = geo.dz[n]; hz[n] = geo.hz[n]; }

		if (n<Nx) { kx[n] = geo.kx[n]; kx2[n] = geo.kx2[n]; }
		if (n<Nz) { kz[n] = geo.kz[n]; kz2[n] = geo.kz2[n]; }

		if (n<=Nx) { ima[n] = geo.ima[n]; ipa[n] = geo.ipa[n]; }
		if (n<=Ny) { jma[n] = geo.jma[n]; jpa[n] = geo.jpa[n]; }
		if (n<=Nz) { kma[n] = geo.kma[n]; kpa[n] = geo.kpa[n]; }
	}
}


void Geometry::InitMesh(double dy_min, const char *path)
{
	for (int i=1; i<=Nx; i++) x[i] = Lx * (i-1.) / (Nx-1);
	for (int k=1; k<=Nz; k++) z[k] = Lz * (k-1.) / (Nz-1);

	if (path) InitMeshY(path);
	else      InitMeshY(dy_min);

	x[0] = 0; // set useless points to 0
	z[0] = 0;
	y[0] = 0;

	for (int i=1; i<Nx; i++) xc[i] = .5 * (x[i] + x[i+1]);
	for (int j=1; j<Ny; j++) yc[j] = .5 * (y[j] + y[j+1]);
	for (int k=1; k<Nz; k++) zc[k] = .5 * (z[k] + z[k+1]);

	xc[0] = x[1]; xc[Nx] = x[Nx];
	yc[0] = y[1]; yc[Ny] = y[Ny];
	zc[0] = z[1]; zc[Nz] = z[Nz];
}

void Geometry::InitInterval()
{
	for (int i=1; i<Nx; i++) dx[i] = x[i+1] - x[i];
	for (int j=1; j<Ny; j++) dy[j] = y[j+1] - y[j];
	for (int k=1; k<Nz; k++) dz[k] = z[k+1] - z[k];

	dx[0] = 2. * (x[1] - xc[0]); // used for interpolation
	dy[0] = 2. * (y[1] - yc[0]);
	dz[0] = 2. * (z[1] - zc[0]);

	dx[Nx] = 2. * (xc[Nx] - x[Nx]);
	dy[Ny] = 2. * (yc[Ny] - y[Ny]);
	dz[Nz] = 2. * (zc[Nz] - z[Nz]);

	for (int i=1; i<=Nx; i++) hx[i] = xc[i] - xc[i-1];
	for (int j=1; j<=Ny; j++) hy[j] = yc[j] - yc[j-1];
	for (int k=1; k<=Nz; k++) hz[k] = zc[k] - zc[k-1];

	hx[0] = 0; // set useless points to 0
	hy[0] = 0;
	hz[0] = 0;
}

void Geometry::InitWaveNumber()
{
	double dx = Lx / (Nx-1);
	double dz = Lz / (Nz-1);

	for (int i=0; i<Nx-1; i++) {
		kx [i] = 2.*PI/Lx * (i - (Nx-1) * (i>=Nxc));
		kx2[i] = 2./dx/dx * (1 - cos(kx[i] * dx));
	}
	for (int k=0; k<Nz-1; k++) {
		kz [k] = 2.*PI/Lz * (k - (Nz-1) * (k>=Nzc));
		kz2[k] = 2./dz/dz * (1 - cos(kz[k] * dz));
	}
}

void Geometry::InitIndices()
{
	for (int i=1; i<Nx; i++) {
		ima[i] = i-1;
		ipa[i] = i+1;
	}
	for (int j=1; j<Ny; j++) {
		jma[j] = j-1;
		jpa[j] = j+1;
	}
	for (int k=1; k<Nz; k++) {
		kma[k] = k-1;
		kpa[k] = k+1;
	}
	// boundary indexing are usually not used, except for calculating taub
	ima[Nx] = Nx-1; ipa[0] = 1;
	jma[Ny] = Ny-1; jpa[0] = 1;
	kma[Nz] = Nz-1; kpa[0] = 1;

	// avoid accessing invalid indices
	ima[0] = ipa[Nx] = 1/INFTSM;
	jma[0] = jpa[Ny] = 1/INFTSM;
	kma[0] = kpa[Nz] = 1/INFTSM;
}




void Geometry::InitMeshY(double dy_min)
/* generate grid points for y-direction
   dy_min >0 for two-sided hyperbolic tangent, <0 one-sided, =0 uniform */
{
	if (fabs(dy_min) < INFTSM) {
		// uniform mesh
		for (int j=1; j<=Ny; j++)
			y[j] = 1 + Ly * ((j-1)/(Ny-1.) - .5);
	}
	else if (fabs(dy_min) <= Ly / (Ny-1.)) {
		// Newton iteration to determine parameter
		double gma = dy_min > 0 ?
			newton_iter( dy_min, get_dymin, Ny, Ly) :
			newton_iter(-dy_min, get_dymin, 2*Ny-1, 2*Ly);

		for (int j=1; j<=Ny; j++)
			y[j] = dy_min > 0 ?
				hyptan(j, gma, Ny, Ly) :
				hyptan(j, gma, 2*Ny-1, 2*Ly); // 1 + Ly * tanh(gma * ((j-1)/(Ny-1) - 1)) / tanh(gma);
	}
	else {
		cout << "Mesh error: dy_min too large !" << endl;
		exit(0);
	}
}

void Geometry::InitMeshY(const char *path)
/* read grid from file at specified path */
{
	FILE *fp;
	char str[1024];

	sprintf(str, "%sCHANNEL.GRD", path);

	fp = fopen(str, "r");

	for (int j=0; j<=Ny; j++)
		sscanf(fgets(str, 1024, fp), "%le", &y[j]);

	fclose(fp);

	// check
	if ( fabs((y[Ny]-y[1]) - Ly) > INFTSM ) {
		cout << endl << "No valid Y mesh provided !" << endl;
		exit(0);
	}
}

void Geometry::AlignBoundaryYc(const Geometry &geo)
// align OFW boundary yc to the FC yc
{
	int jb = (geo.Ny - Ny) / 2;
	yc[0]  = geo.yc[jb];
	yc[Ny] = geo.yc[geo.Ny-jb];
}

void Geometry::WriteMeshY(const char *path) const
{
	char str[1024];
	sprintf(str, "%sCHANNEL.GRD", path);

	FILE *fp = fopen(str, "w");
	for (int j=0; j<=Ny; j++)
		fprintf(fp, "%.18e\n", y[j]);
	fclose(fp);

	// check resolution in y-direction
	double dy_min = 1./INFTSM;
	double dy_max = 0;
	double dy;

	for (int j=1; j<Ny; j++) {
		dy = y[j+1] - y[j];
		if (dy < dy_min) dy_min = dy;
		if (dy > dy_max) dy_max = dy;
	}

	cout << "\ndy_min = " << dy_min
	     << ", dy_max = " << dy_max << endl;
}

void Geometry::WriteMesh(const char *path) const
{
	char str[1024];
	sprintf(str, "%sMESH.txt", path);

	FILE *fp = fopen(str, "w");

	fputs("\nx, xc\n\n", fp);

	for (int i=0; i<=Nx; i++)
		fprintf(fp, "%.18e\t%.18e\n", x[i], xc[i]);

	fputs("\ny, yc\n\n", fp);

	for (int j=0; j<=Ny; j++)
		fprintf(fp, "%.18e\t%.18e\n", y[j], yc[j]);

	fputs("\nz, zc\n\n", fp);

	for (int k=0; k<=Nz; k++)
		fprintf(fp, "%.18e\t%.18e\n", z[k], zc[k]);

	fclose(fp);
}





double hyptan(int i, double gma, int N, double L)
/* generate grid points with hyperbolic tangent distribution */
{
	return 1 + .5*L * tanh(gma * (2.*(i-1)/(N-1) - 1)) / tanh(gma);
	// return .5*L * (1 + tanh(gma * (2.*(i-1)/(N-1) - 1)) / tanh(gma));
}

double get_dymin(double gma, int N, double L)
{
	return hyptan(2, gma, N, L)
		 - hyptan(1, gma, N, L);
}

double newton_iter(double tgt, double (*pfun)(double, int, double), int N, double L)
{
	double F, grad, gma = 1, dgma = .1;
	double F0 = pfun(gma-dgma, N, L) - tgt;

	for (int n=0; n<100; n++) {
		F = pfun(gma, N, L) - tgt;

		grad = (F-F0) / dgma; if (! grad) break;
		dgma = -F / grad;     if (! dgma) break;
		gma += dgma;
		F0 = F;
	}
	return gma;
}






// #define DEBUG

#ifdef DEBUG

int main()
{
	Geometry geo0(4,8,4, 2*PI,2,PI);
	Geometry geo1(4,8,4, 2*PI,1,PI);

	geo0.InitIndices();
	geo1.InitIndices();

	geo0.InitMesh(0);
	Mesh(geo0).WriteMeshY("path_Mesh/test0/");

	geo0.InitMesh(.1);
	Mesh(geo0).WriteMeshY("path_Mesh/test1/");

	geo1.InitMesh(-.02);
	Mesh(geo1).WriteMeshY("path_Mesh/test2/");

	geo1.InitMesh(0, "path_Mesh/test2/");
	Mesh(geo1).WriteMeshY("path_Mesh/test3/");

	return 0;
}

#endif



// // sample compiling commands for test on Mac
// export CPLUS_INCLUDE_PATH=/usr/local/include
// export LIBRARY_PATH=/usr/local/lib
// mkdir -p path_Mesh/{test0,test1,test2,test3}
// g++-9 -lfftw3 -lm -fopenmp Scla.cpp Mesh.cpp -o test_Mesh






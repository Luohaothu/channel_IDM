#include "Basic.h"

using namespace std;


Flow::Flow(const Mesh &ms):
ms(ms),
v_(ms),
s_(ms),
Nx(ms.Nx),
Ny(ms.Ny),
Nz(ms.Nz)
{}


void Flow::InitRand(double energy)
/* initiate velocity (up to real boundaries) with random fluctuations */
{
	Scla &u = v_[1] = 0;
	Scla &v = v_[2] = 0;
	Scla &w = v_[3] = 0;
	Scla &p = s_    = 0;

	srand(time(0));
	for (int j=1; j<=Ny; j++) { // probability distribution: p(x) = ( x<0 ? x+1 : 1-x )
	for (int k=1; k<=Nz; k++) { // 2nd order moment is 2/3, <U^2+V^2+W^2>/2 = energy
	for (int i=1; i<=Nx; i++) { // sample space expands on the whole physical domain
		if (j<Ny && k<Nz) u(i,j,k) = energy * (rand()-rand()) / RAND_MAX;
		if (i<Nx && k<Nz) v(i,j,k) = energy * (rand()-rand()) / RAND_MAX;
		if (i<Nx && j<Ny) w(i,j,k) = energy * (rand()-rand()) / RAND_MAX;
	}}}

	// remove mean velocity normal to every plane
	for (int i=1; i<=Nx; i++) {
		double mean = u.MeanUyz(i);
		for (int j=1; j<Ny; j++)
		for (int k=1; k<Nz; k++)
			u(i,j,k) -= mean;
	}
	for (int j=1; j<=Ny; j++) {
		double mean = v.MeanVxz(j);
		for (int k=1; k<Nz; k++)
		for (int i=1; i<Nx; i++)
			v(i,j,k) -= mean;
	}
	for (int k=1; k<=Nz; k++) {
		double mean = w.MeanWxy(k);
		for (int j=1; j<Ny; j++)
		for (int i=1; i<Nx; i++)
			w(i,j,k) -= mean;
	}
}


/***** file IO operations *****/

Flow& Flow::ReadField(const char *path, int tstep, const char *suffix)
{
	char str[32];
	sprintf(str, "U%s%08i", suffix, tstep); v_[1].FileIO(path, str, 'r');
	sprintf(str, "V%s%08i", suffix, tstep); v_[2].FileIO(path, str, 'r');
	sprintf(str, "W%s%08i", suffix, tstep); v_[3].FileIO(path, str, 'r');
	sprintf(str, "P%s%08i", suffix, tstep);    s_.FileIO(path, str, 'r');
	return *this;
}

void Flow::WriteField(const char *path, int tstep, const char *suffix) const
{
	char str[32];
	sprintf(str, "U%s%08i", suffix, tstep); v_[1].FileIO(path, str, 'w');
	sprintf(str, "V%s%08i", suffix, tstep); v_[2].FileIO(path, str, 'w');
	sprintf(str, "W%s%08i", suffix, tstep); v_[3].FileIO(path, str, 'w');
	sprintf(str, "P%s%08i", suffix, tstep);    s_.FileIO(path, str, 'w');
}

void Flow::WriteTecplot(const char *path, int tstep, double time) const
/* write velocity and pressure fields to ascii files readable by tecplot */
{
	FILE *fp;
	char str[1024], filename[1024];

	const Scla &p = s_;
	Scla u(ms);
	Scla v(ms);
	Scla w(ms);

	v_[1].Ugrid2CellCenter(u);
	v_[2].Vgrid2CellCenter(v);
	v_[3].Wgrid2CellCenter(w);

	sprintf(filename, "%sFIELD%08i.dat", path, tstep);

	fp = fopen(filename, "w");

	fputs("Title = \"3D instantaneous field\"\n", fp);
	fputs("variables = \"x\", \"y\", \"z\", \"u\", \"v\", \"w\", \"p\"\n", fp);
	fprintf(fp, "zone t = \"%f\", i = %i, j = %i, k = %i\n", time, Nx+1, Nz+1, Ny+1);

	for (int j=0; j<=Ny; j++) {
	for (int k=0; k<=Nz; k++) {
	for (int i=0; i<=Nx; i++) {
		
		fprintf(fp,
			"%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
			ms.xc(i),
			ms.yc(j),
			ms.zc(k),
			u(i,j,k),
			v(i,j,k),
			w(i,j,k),
			p(i,j,k) );
	}}}

	fclose(fp);

	sprintf(str, "preplot %s", filename); system(str);
	sprintf(str, "rm %s",      filename); system(str);
}












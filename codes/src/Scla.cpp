#include "Basic.h"

using namespace std;

typedef fftw_complex fcmplx;



Scla::Scla(const Mesh &ms):
ms(ms),
Nx(ms.Nx),
Ny(ms.Ny),
Nz(ms.Nz)
{
	q_ = new double[(Nx+1) * (Ny+1) * (Nz+1)];
	memset(q_, 0, (Nx+1) * (Ny+1) * (Nz+1) * sizeof(double));

	fft_temp = new double*[Ny+1];
	for (int j=0; j<=Ny; j++)
		fft_temp[j] = new double[max(ms.Nzr,ms.Nxr)];

	frc_x = new fftw_plan[Ny+1];
	fcr_x = new fftw_plan[Ny+1];
	frc_z = new fftw_plan[Ny+1];
	fcr_z = new fftw_plan[Ny+1];
	frc_xz= new fftw_plan[Ny+1];
	fcr_xz= new fftw_plan[Ny+1];
	frR_x = new fftw_plan[(Ny+1)*(Nz+1)];
	fRr_x = new fftw_plan[(Ny+1)*(Nz+1)];

	for (int j=0; j<=Ny; j++) {

		double *r = fft_temp[j];
		fcmplx *c = (fcmplx*) r;

		frc_x[j] = fftw_plan_dft_r2c_1d(Nx-1, r, c, FFTW_MEASURE);
		fcr_x[j] = fftw_plan_dft_c2r_1d(Nx-1, c, r, FFTW_MEASURE);

		frc_z[j] = fftw_plan_dft_r2c_1d(Nz-1, r, c, FFTW_MEASURE);
		fcr_z[j] = fftw_plan_dft_c2r_1d(Nz-1, c, r, FFTW_MEASURE);

		r = &q_[ms.idfxz(0,j,0)];
		c = (fcmplx*) r;

		frc_xz[j] = fftw_plan_dft_r2c_2d(Nz-1, Nx-1, r, c, FFTW_MEASURE);
		fcr_xz[j] = fftw_plan_dft_c2r_2d(Nz-1, Nx-1, c, r, FFTW_MEASURE);

		for (int k=0; k<=Nz; k++) {
			r = &q_[ms.idx(1,j,k)];

			frR_x[j*(Nz+1)+k] = fftw_plan_r2r_1d(Nx-1, r, r, FFTW_REDFT10, FFTW_MEASURE);
			fRr_x[j*(Nz+1)+k] = fftw_plan_r2r_1d(Nx-1, r, r, FFTW_REDFT01, FFTW_MEASURE);
		}
	}
}

Scla::~Scla()
{
	for (int j=0; j<=Ny; j++) { // note: must not destory fft plan that is not created
		fftw_destroy_plan(frc_x[j]);
		fftw_destroy_plan(fcr_x[j]);
		fftw_destroy_plan(frc_z[j]);
		fftw_destroy_plan(fcr_z[j]);
		fftw_destroy_plan(frc_xz[j]);
		fftw_destroy_plan(fcr_xz[j]);
		for (int k=0; k<=Nz; k++) {
			fftw_destroy_plan(frR_x[j*(Nz+1)+k]);
			fftw_destroy_plan(fRr_x[j*(Nz+1)+k]);
		}
		delete[] fft_temp[j];
	}
	delete[] frc_x;
	delete[] fcr_x;
	delete[] frc_z;
	delete[] fcr_z;
	delete[] frc_xz;
	delete[] fcr_xz;
	delete[] frR_x;
	delete[] fRr_x;
	delete[] fft_temp;
	delete[] q_;
}

Scla::Scla(const Scla &src):
Scla(src.ms)
{
	Set(src);
}

Scla& Scla::operator=(const Scla &src)
{
	if (&src == this) return *this;
	if (src.ms.Nx != ms.Nx ||
		src.ms.Ny != ms.Ny ||
		src.ms.Nz != ms.Nz)
	{
		cout << "Scalar sizes do not match !" << endl;
		exit(0);
	}
	cout << "Warning: if you just want to assign values, use Set() instead." << endl;
	return Set(src);
}


/***** fft *****/

void Scla::fftxz(int j0, int jn)
{
	#pragma omp for
	for (int j=j0; j<=(jn>0?jn:Ny+jn); j++) {
		for (int k=1; k<Nz; k++)
		for (int i=1; i<Nx; i++)
			q_[ms.idfxz(i-1,j,k-1)] = q_[ms.idx(i,j,k)];
		fftw_execute(frc_xz[j]);
	}
}
void Scla::ifftxz(int j0, int jn)
{
	#pragma omp for
	for (int j=j0; j<=(jn>0?jn:Ny+jn); j++) {
		fftw_execute(fcr_xz[j]);
		for (int k=Nz-1; k>=1; k--)
		for (int i=Nx-1; i>=1; i--)
			q_[ms.idx(i,j,k)] = q_[ms.idfxz(i-1,j,k-1)] / (Nx-1)/(Nz-1);
	}
}

void Scla::fftx(int j0, int jn)
{
	#pragma omp for
	for (int j=j0; j<=(jn>0?jn:Ny+jn); j++) { double *temp = fft_temp[j];
	for (int k=0; k<=Nz; k++) {
		for (int i=1; i<Nx; i++)
			temp[i-1] = q_[ms.idx(i,j,k)];
		fftw_execute(frc_x[j]);
		for (int i=0; i<ms.Nxr; i++)
			q_[ms.idx(i,j,k)] = temp[i];
	}}
}
void Scla::ifftx(int j0, int jn)
{
	#pragma omp for
	for (int j=j0; j<=(jn>0?jn:Ny+jn); j++) { double *temp = fft_temp[j];
	for (int k=0; k<=Nz; k++) {
		for (int i=0; i<ms.Nxr; i++)
			temp[i] = q_[ms.idx(i,j,k)];
		fftw_execute(fcr_x[j]);
		for (int i=1; i<Nx; i++)
			q_[ms.idx(i,j,k)] = temp[i-1] / (Nx-1);
	}}
}

void Scla::fftz(int j0, int jn)
{
	#pragma omp for
	for (int j=j0; j<=(jn>0?jn:Ny+jn); j++) { double *temp = fft_temp[j];
	for (int i=0; i<=Nx; i++) {
		for (int k=1; k<Nz; k++)
			temp[k-1] = q_[ms.idx(i,j,k)];
		fftw_execute(frc_z[j]);
		for (int k=0; k<ms.Nzr; k++)
			q_[ms.idx(i,j,k)] = temp[k];
	}}
}
void Scla::ifftz(int j0, int jn)
{
	#pragma omp for
	for (int j=j0; j<=(jn>0?jn:Ny+jn); j++) { double *temp = fft_temp[j];
	for (int i=0; i<=Nx; i++) {
		for (int k=0; k<ms.Nzr; k++)
			temp[k] = q_[ms.idx(i,j,k)];
		fftw_execute(fcr_z[j]);
		for (int k=1; k<Nz; k++)
			q_[ms.idx(i,j,k)] = temp[k-1] / (Nz-1);
	}}
}

void Scla::dctxz(int j0, int jn)
{
	this->dctx(j0, jn);
	this->fftz(j0, jn);
}

void Scla::idctxz(int j0, int jn)
{
	this->ifftz(j0, jn);
	this->idctx(j0, jn);
}

void Scla::dctx(int j0, int jn)
{
	#pragma omp for
	for (int j=j0; j<=(jn>0?jn:Ny+jn); j++) {
	for (int k=0; k<=Nz; k++) {
		fftw_execute(frR_x[j*(Nz+1)+k]);
		for (int i=1; i<Nx; i++)
			q_[ms.idx(i,j,k)] /= 2;
	}}
}
void Scla::idctx(int j0, int jn)
{
	#pragma omp for
	for (int j=j0; j<=(jn>0?jn:Ny+jn); j++) {
	for (int k=0; k<=Nz; k++) {
		fftw_execute(fRr_x[j*(Nz+1)+k]);
		for (int i=1; i<Nx; i++)
			q_[ms.idx(i,j,k)] /= Nx-1;
	}}
}


/***** convinent operations for whole arrays *****/

void Scla::TravLyr(double a, int j, void (*pfun)(double &b, double a))
{
	for (int k=0; k<=Nz; k++)
	for (int i=0; i<=Nx; i++)
		pfun(q_[ms.idx(i,j,k)], a);
}
void Scla::TravLyr(const double *src, int j, void (*pfun)(double &b, double a))
{
	for (int k=0; k<=Nz; k++)
	for (int i=0; i<=Nx; i++)
		pfun(q_[ms.idx(i,j,k)], src[ms.idx(i,0,k)]);
}
void Scla::TravBlk(double a, void (*pfun)(double&, double))
{
	#pragma omp parallel for
	for (int j=0; j<=Ny; j++)
	for (int k=0; k<=Nz; k++)
	for (int i=0; i<=Nx; i++)
		pfun(q_[ms.idx(i,j,k)], a);
}
void Scla::TravBlk(const double *src, void (*pfun)(double&, double))
{
	#pragma omp parallel for
	for (int j=0; j<=Ny; j++) {
	for (int k=0; k<=Nz; k++) {
	for (int i=0; i<=Nx; i++) {
		int id = ms.idx(i,j,k);
		pfun(q_[id], src[id]);
	}}}
}


/***** arithmetic mean *****/

double Scla::meanxz(int j) const
{
	double sum = 0;
	for (int k=1; k<Nz; k++)
	for (int i=1; i<Nx; i++)
		sum += q_[ms.idx(i,j,k)];
	return sum / ((Nx-1)*(Nz-1));
}
double Scla::meanz(int i, int j) const
{
	double sum = 0;
	for (int k=1; k<Nz; k++)
		sum += q_[ms.idx(i,j,k)];
	return sum / (Nz-1);
}


/***** weighed average *****/

double Scla::MeanUx(int j, int k) const
{
	double sum = 0;
	for (int i=2; i<Nx; i++)
		sum += q_[ms.idx(i,j,k)] * ms.hx(i);
	sum += q_[ms.idx(1, j,k)] * (ms.xc(1) - ms.x(1));
	sum += q_[ms.idx(Nx,j,k)] * (ms.x(Nx) - ms.xc(Nx-1));
	return sum / ms.Lx;
}
double Scla::MeanVy(int i, int k) const
{
	double sum = 0;
	for (int j=2; j<Ny; j++)
		sum += q_[ms.idx(i,j,k)] * ms.hy(j);
	sum += q_[ms.idx(i,1, k)] * (ms.yc(1) - ms.y(1));
	sum += q_[ms.idx(i,Ny,k)] * (ms.y(Ny) - ms.yc(Ny-1));
	return sum / ms.Ly;
}
double Scla::MeanWz(int i, int j) const
{
	double sum = 0;
	for (int k=2; k<Nz; k++)
		sum += q_[ms.idx(i,j,k)] * ms.hz(k);
	sum += q_[ms.idx(i,j,1 )] * (ms.zc(1) - ms.z(1));
	sum += q_[ms.idx(i,j,Nz)] * (ms.z(Nz) - ms.zc(Nz-1));
	return sum / ms.Lz;
}

double Scla::MeanAy(int i, int k) const
{
	double sum = 0;
	for (int j=1; j<Ny; j++)
		sum += q_[ms.idx(i,j,k)] * ms.dy(j);
	return sum / ms.Ly;
}

double Scla::MeanAz(int i, int j) const
{
	double sum = 0;
	for (int k=1; k<Nz; k++)
		sum += q_[ms.idx(i,j,k)] * ms.dz(k);
	return sum / ms.Lz;
}

double Scla::MeanUyz(int i) const
{
	double sum = 0;
	for (int j=1; j<Ny; j++)
	for (int k=1; k<Nz; k++)
		sum += q_[ms.idx(i,j,k)] * ms.dy(j) * ms.dz(k);
	return sum / (ms.Ly * ms.Lz);
}
double Scla::MeanVxz(int j) const
{
	double sum = 0;
	for (int k=1; k<Nz; k++)
	for (int i=1; i<Nx; i++)
		sum += q_[ms.idx(i,j,k)] * ms.dx(i) * ms.dz(k);
	return sum / (ms.Lx * ms.Lz);

}
double Scla::MeanWxy(int k) const
{
	double sum = 0;
	for (int j=1; j<Ny; j++)
	for (int i=1; i<Nx; i++)
		sum += q_[ms.idx(i,j,k)] * ms.dx(i) * ms.dy(j);
	return sum / (ms.Lx * ms.Ly);
}

double Scla::MeanU() const
{
	double sum = 0;
	#pragma omp parallel for reduction(+: sum)
	for (int j=1; j<Ny; j++)
	for (int k=1; k<Nz; k++)
		sum += MeanUx(j,k) * ms.dy(j) * ms.dz(k);
	return sum / (ms.Ly * ms.Lz);
}
double Scla::MeanV() const
{
	double sum = 0;
	#pragma omp parallel for reduction(+: sum)
	for (int k=1; k<Nz; k++)
	for (int i=1; i<Nx; i++)
		sum += MeanVy(i,k) * ms.dx(i) * ms.dz(k);
	return sum / (ms.Lx * ms.Lz);
}
double Scla::MeanW() const
{
	double sum = 0;
	#pragma omp parallel for reduction(+: sum)
	for (int j=1; j<Ny; j++)
	for (int i=1; i<Nx; i++)
		sum += MeanWz(i,j) * ms.dx(i) * ms.dy(j);
	return sum / (ms.Lx * ms.Ly);
}
double Scla::MeanA() const
{
	double sum = 0;
	#pragma omp parallel for reduction(+: sum)
	for (int j=1; j<Ny; j++)
		sum += MeanVxz(j) * ms.dy(j);
	return sum / ms.Ly;
}


/***** interpolate from faces to cell-centers *****/

// way to deal with real-to-virtual booundary interpolation:
// linear extrapolation, without special consideration for periodic scenario

Scla& Scla::Ugrid2CellCenter(const Scla &src)
{
	double c1 = .5 * ms.dx(0) / ms.dx(1);
	double c2 = .5 * ms.dx(Nx)/ ms.dx(Nx-1);

	#pragma omp for collapse(2)
	for (int j=0; j<=Ny; j++) {
	for (int k=0; k<=Nz; k++) {
		for (int i=1; i<Nx; i++) {
			q_[ms.idx(i,j,k)] = .5 * (src(i,j,k) + src(ms.ipa(i),j,k));
		}
		q_[ms.idx(0, j,k)] = (1+c1) * src(1, j,k) - c1 * src(2,j,k);
		q_[ms.idx(Nx,j,k)] = (1+c2) * src(Nx,j,k) - c2 * src(Nx-1,j,k);
	}}
	return *this;
}
Scla& Scla::Vgrid2CellCenter(const Scla &src)
{
	double c1 = .5 * ms.dy(0) / ms.dy(1);
	double c2 = .5 * ms.dy(Ny)/ ms.dy(Ny-1);

	#pragma omp for collapse(2)
	for (int k=0; k<=Nz; k++) {
	for (int i=0; i<=Nx; i++) {
		for (int j=1; j<Ny; j++) {
			q_[ms.idx(i,j,k)] = .5 * (src(i,j,k) + src(i,ms.jpa(j),k));
		}
		q_[ms.idx(i,0, k)] = (1+c1) * src(i,1, k) - c1 * src(i,2,k);
		q_[ms.idx(i,Ny,k)] = (1+c2) * src(i,Ny,k) - c2 * src(i,Ny-1,k);
	}}
	return *this;
}
Scla& Scla::Wgrid2CellCenter(const Scla &src)
{
	double c1 = .5 * ms.dz(0) / ms.dz(1);
	double c2 = .5 * ms.dz(Nz)/ ms.dz(Nz-1);

	#pragma omp for collapse(2)
	for (int j=0; j<=Ny; j++) {
	for (int i=0; i<=Nx; i++) {
		for (int k=1; k<Nz; k++) {
			q_[ms.idx(i,j,k)] = .5 * (src(i,j,k) + src(i,j,ms.kpa(k)));
		}
		q_[ms.idx(i,j,0 )] = (1+c1) * src(i,j,1 ) - c1 * src(i,j,2);
		q_[ms.idx(i,j,Nz)] = (1+c2) * src(i,j,Nz) - c2 * src(i,j,Nz-1);
	}}
	return *this;
}


/***** interpolate from cell-centers to edges *****/

void Scla::CellCenter2Edges(Scla &dst1, Scla &dst2, Scla &dst3) const
{
	#pragma omp for
	for (int j=0; j<=Ny; j++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hyc=ms.hy(j);
	for (int k=0; k<=Nz; k++) { double dzm,dzp,dzc=ms.dz(k,dzm,dzp), hzc=ms.hz(k);
	for (int i=0; i<=Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxc=ms.hx(i);

		int id =              ms.idx(i,j,k);
		int im, jm, km;       ms.imx(i,j,k,im,jm,km);
		int imjm, jmkm, imkm; ms.mmx(i,j,k,imjm,jmkm,imkm);

		if (j>0 && k>0) dst1[id] = .25/hyc/hzc * (
			dym * dzm * q_[id] + dyc * dzm * q_[jm] +
			dym * dzc * q_[km] + dyc * dzc * q_[jmkm] );

		if (i>0 && k>0) dst2[id] = .25/hxc/hzc * (
			dxm * dzm * q_[id] + dxc * dzm * q_[im] +
			dxm * dzc * q_[km] + dxc * dzc * q_[imkm] );

		if (j>0 && i>0) dst3[id] = .25/hxc/hyc * (
			dxm * dym * q_[id] + dxc * dym * q_[im] +
			dxm * dyc * q_[jm] + dxc * dyc * q_[imjm] );
	}}}
}

/***** differentiation operators *****/

vector<double>  Scla::Gradient(int i, int j, int k) const
/* compute Gradient of a cell-centered scalar field to corresponding U,V,W grids */

{
	vector<double> grad(3,0.);

	int id =        ms.idx(i,j,k);
	int im, jm, km; ms.imx(i,j,k,im,jm,km);

	grad[0] = (q_[id] - q_[im]) / ms.hx(i); // [1,Nx], [0,Ny], [0,Nz]
	grad[1] = (q_[id] - q_[jm]) / ms.hy(j); // [0,Nx], [1,NY], [0,Nz]
	grad[2] = (q_[id] - q_[km]) / ms.hz(k); // [0,Nx], [0,NY], [1,Nz]

	return grad;
}


/***** IO *****/

void Scla::FileIO(const char *path, const char *name, char mode) const
/* read & write field from & to binary files */
{
	FILE *fp;
	char str[1024];

	sprintf(str, "%s%s.bin", path, name);

	fp = fopen(str, mode=='w' ? "wb" : "rb");

	// write domain information at the beginning
	int n1 = Nx+1;
	int n2 = Nz+1;
	int n3 = Ny+1;
	if (mode == 'w') {
		fwrite(&n1, sizeof(int), 1, fp);
		fwrite(&n2, sizeof(int), 1, fp);
		fwrite(&n3, sizeof(int), 1, fp);
	}
	// data begin after the info section
	fseek(fp, sizeof(double) * n1*n2, SEEK_SET);
	if (mode == 'w') fwrite(q_, sizeof(double) * n1*n2, n3, fp);
	if (mode == 'r') fread (q_, sizeof(double) * n1*n2, n3, fp);

	fclose(fp);
}

void Scla::debug_AsciiOutput(const char *path, const char *name, int j1, int j2) const
/* write the fields in ascii files for check */
{
	FILE *fp;
	char str[1024];

	sprintf(str, "%s%s.txt", path, name);

	fp = fopen(str, "w");
	for (int j=j1;j<j2;  j++) { fputc('\n', fp);
	for (int k=0; k<=Nz; k++) { fputc('\n', fp);
	for (int i=0; i<=Nx; i++) {
		fprintf(fp, "%.6f\t", q_[ms.idx(i,j,k)]);
	}}}
	fclose(fp);
}





// #define DEBUG

#ifdef DEBUG

int main()
{
	Geometry geo(5,3,4, 2*PI,2,PI);
	Mesh ms(geo);
	Scla q(ms);

	geo.InitIndices();
	geo.InitMesh(.2);

	for (int j=1; j<Ny; j++) {
	for (int k=1; k<Nz; k++) {
	for (int i=1; i<Nx; i++) {
		q(i,j,k) = i*j + j*k + k*i;
	}}}

	q.debug_AsciiOutput("", "0_original", 0, Ny+1);
	q.fftz();
	q.debug_AsciiOutput("", "1_fftz", 0, Ny+1);
	q.ifftz();
	q.fftxz();
	q.debug_AsciiOutput("", "2_fftxz", 0, Ny+1);
	q.ifftxz();
	q.debug_AsciiOutput("", "3_back", 0, Ny+1);

	return 0;
}

#endif



// // sample compiling commands for test on Mac
// export CPLUS_INCLUDE_PATH=/usr/local/include
// export LIBRARY_PATH=/usr/local/lib
// g++-9 -lfftw3 -lm -fopenmp Scla.cpp Mesh.cpp -o test_Scla



// // sample verifying codes in Python
// import numpy as np
// from numpy.fft import fft, fft2

// Nx, Ny, Nz = 5, 3, 4

// q = np.zeros(Ny+1, Nz+1, Nx+1)
// for j in range(1, Ny):
// 	for k in range(1, Nz):
// 		for i in range(1, Nx):
// 			q[j,k,i] = i*j + j*k + k*i

// for j in range(0, Ny+1):
// 	print(fft(q[j,1:-1], axis=-2))

// for j in range(0, Ny+1):
// 	print(fft2(q[j,1:-1,1:-1], axis=(-1,-2)))




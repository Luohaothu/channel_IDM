#include "Filter.h"
#include "Interp.h"

using namespace std;


static double BoxFilter(double xm, double xp, const double *xs, const double *f, int i0, int in)
// box filter in range [xm,xp] of function f(x) defined on [xs_i0,xs_in]
{
	// handle periodic scenario, assume f[i0] = f[in]
	if (xm > xp) {
		double x1 = xs[in] - xm;
		double x2 = xp - xs[i0];
		return 1. / (x1 + x2) * (
			x1 * BoxFilter(xm, xs[in], xs, f, i0, in)
		  + x2 * BoxFilter(xs[i0], xp, xs, f, i0, in) );
	}

	// binary search for the interpolation range
	int im = Interp::BiSearch(xm, xs, i0, in);     // im \in [i0, in-1]
	int ip = Interp::BiSearch(xp, xs, i0, in) + 1; // ip \in [i0+1, in]

	double sum = 0;

	double x1 = xm;
	double x2 = xs[im+1];

	double f2 = f[im+1];
	double f1 = (f[im] * (x2-xm) + f2 * (xm-xs[im])) / (x2 - xs[im]);

	int i = im;

	while (xs[++i] < xp) {

		x2 = xs[i];
		f2 = f[i];

		sum += .5 * (f1+f2) * (x2-x1);

		x1 = x2;
		f1 = f2;
	}

	x2 = xp;
	f2 = (f1 * (xs[ip]-xp) + f[ip] * (xp-x1)) / (xs[ip] - x1);
	sum += .5 * (f1+f2) * (x2-x1);

	return sum / (xp-xm);
}

static double BoxFilter2(double xm, double xp, double ym, double yp,
	const double *xs, const double *ys, const double **f, int i0, int in, int j0, int jn)
// box filter in range [xm,xp]*[ym,yp] of function f(x,y) defined on [xs_i0,xs_in]*[ys_j0,ys_jn]
{
	// binary search for the interpolation range
	int i1 = Interp::BiSearch(xm, xs, i0, in);     // i1 \in [i0, in-1]
	int i2 = Interp::BiSearch(xp, xs, i0, in) + 1; // i2 \in [i0+1, in]
	int j1 = Interp::BiSearch(ym, ys, j0, jn);     // j1 \in [j0, jn-1]
	int j2 = Interp::BiSearch(yp, ys, j0, jn) + 1; // j2 \in [j0+1, jn]

	// handle periodic scenario, assume f[:,i0]=f[:,in], f[j0,:]=f[jn,:]
	if (ym > yp) {
		double y1 = ys[jn] - ym;
		double y2 = yp - ys[j0];
		return 1. / (y1 + y2) * (
			y1 * BoxFilter2(xm,xp,ym,ys[jn], xs,ys,f, i0,in,j1,jn) +
			y2 * BoxFilter2(xm,xp,ys[j0],yp, xs,ys,f, i0,in,j0,j2) );
	}
	else if (xm > xp) {
		double x1 = xs[in] - xm;
		double x2 = xp - xs[i0];
		return 1. / (x1 + x2) * (
			x1 * BoxFilter2(xm,xs[in],ym,yp, xs,ys,f, i1,in,j0,jn) +
			x2 * BoxFilter2(xs[i0],xp,ym,yp, xs,ys,f, i0,i2,j0,jn) );
	}

	// divide problem into smaller ones
	if (j2-j1 > 1) {
		int j = (j1+j2) / 2;
		double y = ys[j];
		double y1 = y - ym;
		double y2 = yp - y;
		return 1. / (y1 + y2) * (
			y1 * BoxFilter2(xm,xp,ym,y, xs,ys,f, i0,in,j0,j) +
			y2 * BoxFilter2(xm,xp,y,yp, xs,ys,f, i0,in,j,jn) );
	}
	else if (i2-i1 > 1) {
		int i = (i1+i2) / 2;
		double x = xs[i];
		double x1 = x - xm;
		double x2 = xp - x;
		return 1. / (x1 + x2) * (
			x1 * BoxFilter2(xm,x,ym,yp, xs,ys,f, i0,i,j0,jn) +
			x2 * BoxFilter2(x,xp,ym,yp, xs,ys,f, i,in,j0,jn) );
	}

	// tackle only one cell
	// *********************** //
	//       i1            i2
	//       x1  xm    xp  x2
	// j1 y1  -------------
	//       |             |
	//    ym |  fmm---fmp  |
	//       |   |     |   |
	//    yp |  fpm---fpp  |
	//       |             |
	// j2 y2  -------------
	// *********************** //

	double x1 = xs[i1], x2 = xs[i2];
	double y1 = ys[j1], y2 = ys[j2];

	double fmm = 1. / (x2-x1) / (y2-y1) * (
		f[j1][i1] * (x2-xm) * (y2-ym) +
		f[j1][i2] * (xm-x1) * (y2-ym) +
		f[j2][i1] * (x2-xm) * (ym-y1) +
		f[j2][i2] * (xm-x1) * (ym-y1) );

	double fmp = 1. / (x2-x1) / (y2-y1) * (
		f[j1][i1] * (x2-xp) * (y2-ym) +
		f[j1][i2] * (xp-x1) * (y2-ym) +
		f[j2][i1] * (x2-xp) * (ym-y1) +
		f[j2][i2] * (xp-x1) * (ym-y1) );

	double fpm = 1. / (x2-x1) / (y2-y1) * (
		f[j1][i1] * (x2-xm) * (y2-yp) +
		f[j1][i2] * (xm-x1) * (y2-yp) +
		f[j2][i1] * (x2-xm) * (yp-y1) +
		f[j2][i2] * (xm-x1) * (yp-y1) );

	double fpp = 1. / (x2-x1) / (y2-y1) * (
		f[j1][i1] * (x2-xp) * (y2-yp) +
		f[j1][i2] * (xp-x1) * (y2-yp) +
		f[j2][i1] * (x2-xp) * (yp-y1) +
		f[j2][i2] * (xp-x1) * (yp-y1) );

	// automatically degenerate to interpolation if xm=xp and ym=yp
	return .25 * (fmm+fmp+fpm+fpp);
}

static double BoxFilter3(
	double xm, double xp, double ym, double yp, double zm, double zp,
	const double *xs, const double *ys, const double *zs, const Scla &q,
	int i0, int in, int j0, int jn, int k0, int kn)
// box filter in range [xm,xp]*[ym,yp]*[zm,zp] of function f(x,y,z)
// defined on [xs_i0,xs_in]*[ys_j0,ys_jn]*[zs_k0,zs_kn]
{
	// binary search for the interpolation range
	int i1 = Interp::BiSearch(xm, xs, i0, in);     // i1 \in [i0, in-1]
	int i2 = Interp::BiSearch(xp, xs, i0, in) + 1; // i2 \in [i0+1, in]
	int j1 = Interp::BiSearch(ym, ys, j0, jn);     // j1 \in [j0, jn-1]
	int j2 = Interp::BiSearch(yp, ys, j0, jn) + 1; // j2 \in [j0+1, jn]
	int k1 = Interp::BiSearch(zm, zs, k0, kn);     // k1 \in [k0, kn-1]
	int k2 = Interp::BiSearch(zp, zs, k0, kn) + 1; // k2 \in [k0+1, kn]

	// handle periodic scenario, assume f[:,:,i0]=f[:,:,in], f[:,j0,:]=f[:,jn,:], f[k0,:,:]=f[kn,:,:]
	if (zm > zp) {
		double z1 = zs[kn] - zm;
		double z2 = zp - zs[k0];
		return 1. / (z1 + z2) * (
			z1 * BoxFilter3(xm,xp,ym,yp,zm,zs[kn], xs,ys,zs,q, i0,in,j0,jn,k1,kn) +
			z2 * BoxFilter3(xm,xp,ym,yp,zs[k0],zp, xs,ys,zs,q, i0,in,j0,jn,k0,k2) );
	}
	else if (ym > yp) {
		double y1 = ys[jn] - ym;
		double y2 = yp - ys[j0];
		return 1. / (y1 + y2) * (
			y1 * BoxFilter3(xm,xp,ym,ys[jn],zm,zp, xs,ys,zs,q, i0,in,j1,jn,k0,kn) +
			y2 * BoxFilter3(xm,xp,ys[j0],yp,zm,zp, xs,ys,zs,q, i0,in,j0,j2,k0,kn) );
	}
	else if (xm > xp) {
		double x1 = xs[in] - xm;
		double x2 = xp - xs[i0];
		return 1. / (x1 + x2) * (
			x1 * BoxFilter3(xm,xs[in],ym,yp,zm,zp, xs,ys,zs,q, i1,in,j0,jn,k0,kn) +
			x2 * BoxFilter3(xs[i0],xp,ym,yp,zm,zp, xs,ys,zs,q, i0,i2,j0,jn,k0,kn) );
	}

	int recursdep = (i2-i1) * (j2-j1) * (k2-k1); // lower value for deeper recursion
	int threshold = 32; // avoid over-parallism

	int i, j, k;
	double x, y, z, f1, f2;

	// divide and conquer
	if (j2-j1 > 1) {
		y = ys[j = (j1+j2) / 2];
		// #pragma omp task shared(f1, xm,xp,ym,y,zm,zp, xs,ys,zs,q, i0,in,j0,j,k0,kn) if (recursdep > threshold)
		f1 = BoxFilter3(xm,xp,ym,y,zm,zp, xs,ys,zs,q, i0,in,j0,j,k0,kn); // only line that is tasked by omp
		f2 = BoxFilter3(xm,xp,y,yp,zm,zp, xs,ys,zs,q, i0,in,j,jn,k0,kn);
		// #pragma omp taskwait
		return 1. / (yp-ym) * ((y-ym) * f1 + (yp-y) * f2);
	}
	else if (k2-k1 > 1) {
		z = zs[k = (k1+k2) / 2];
		// #pragma omp task shared(f1, xm,xp,ym,yp,zm,z, xs,ys,zs,q, i0,in,j0,jn,k0,k) if (recursdep > threshold)
		f1 = BoxFilter3(xm,xp,ym,yp,zm,z, xs,ys,zs,q, i0,in,j0,jn,k0,k); // only line that is tasked by omp
		f2 = BoxFilter3(xm,xp,ym,yp,z,zp, xs,ys,zs,q, i0,in,j0,jn,k,kn);
		// #pragma omp taskwait
		return 1. / (zp-zm) * ((z-zm) * f1 + (zp-z) * f2);
	}
	else if (i2-i1 > 1) {
		x = xs[i = (i1+i2) / 2];
		// #pragma omp task shared(f1, xm,x,ym,yp,zm,zp, xs,ys,zs,q, i0,i,j0,jn,k0,kn) if (recursdep > threshold)
		f1 = BoxFilter3(xm,x,ym,yp,zm,zp, xs,ys,zs,q, i0,i,j0,jn,k0,kn); // only line that is tasked by omp
		f2 = BoxFilter3(x,xp,ym,yp,zm,zp, xs,ys,zs,q, i,in,j0,jn,k0,kn);
		// #pragma omp taskwait
		return 1. / (xp-xm) * ((x-xm) * f1 + (xp-x) * f2);
	}

	double x1 = xs[i1], x2 = xs[i2];
	double y1 = ys[j1], y2 = ys[j2];
	double z1 = zs[k1], z2 = zs[k2];

	double fmmm = Interp::InterpCell(q, xm,x1,x2,i1,i2, ym,y1,y2,j1,j2, zm,z1,z2,k1,k2);
	double fmmp = Interp::InterpCell(q, xp,x1,x2,i1,i2, ym,y1,y2,j1,j2, zm,z1,z2,k1,k2);
	double fmpm = Interp::InterpCell(q, xm,x1,x2,i1,i2, yp,y1,y2,j1,j2, zm,z1,z2,k1,k2);
	double fmpp = Interp::InterpCell(q, xp,x1,x2,i1,i2, yp,y1,y2,j1,j2, zm,z1,z2,k1,k2);
	double fpmm = Interp::InterpCell(q, xm,x1,x2,i1,i2, ym,y1,y2,j1,j2, zp,z1,z2,k1,k2);
	double fpmp = Interp::InterpCell(q, xp,x1,x2,i1,i2, ym,y1,y2,j1,j2, zp,z1,z2,k1,k2);
	double fppm = Interp::InterpCell(q, xm,x1,x2,i1,i2, yp,y1,y2,j1,j2, zp,z1,z2,k1,k2);
	double fppp = Interp::InterpCell(q, xp,x1,x2,i1,i2, yp,y1,y2,j1,j2, zp,z1,z2,k1,k2);

	return .125 * (fmmm+fmmp+fmpm+fmpp+fpmm+fpmp+fppm+fppp);
}



double Filter::FilterNode(
	double x,  double y,  double z,
	double dx, double dy, double dz, const Scla &q, int stgtyp)
// Box filter at physical point x_i in neighbourhood dx_i (not necessarily = filter size)
// In directions where dx_i = 0, automatically degenerate to linear interpolation
{
	const Mesh &ms = q.ms;

	// configure the staggered coordinate to be interpolated
	const double *corx = stgtyp==1 ? ms.x() : ms.xc();
	const double *cory = stgtyp==2 ? ms.y() : ms.yc();
	const double *corz = stgtyp==3 ? ms.z() : ms.zc();

	// configure the effective domain (0~n for non-periodic xc, 1~n otherwise)
	int in = ms.Nx, i0 = (stgtyp==1 || corx[in]-corx[0]>ms.Lx+INFTSM);
	int jn = ms.Ny, j0 = (stgtyp==2 || cory[jn]-cory[0]>ms.Ly+INFTSM);
	int kn = ms.Nz, k0 = (stgtyp==3 || corz[kn]-corz[0]>ms.Lz+INFTSM);

	// shift target position into meaningful range
	double xm = Interp::Shift(x-.5*dx, corx[i0], corx[in]);
	double ym = Interp::Shift(y-.5*dy, cory[j0], cory[jn]);
	double zm = Interp::Shift(z-.5*dz, corz[k0], corz[kn]);

	double xp = Interp::Shift(x+.5*dx, corx[i0], corx[in]);
	double yp = Interp::Shift(y+.5*dy, cory[j0], cory[jn]);
	double zp = Interp::Shift(z+.5*dz, corz[k0], corz[kn]);

	double ans;

	// #pragma omp parallel
	// #pragma omp single
	ans = BoxFilter3(xm,xp,ym,yp,zm,zp, corx,cory,corz,q, i0,in,j0,jn,k0,kn);
	
	return ans;
}

double Filter::TestFilter(int i, int j, int k, const Scla &q)
{
	const Mesh &ms = q.ms;

	int id =              ms.idx(i,j,k);
	int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
	int im, jm, km;       ms.imx(i,j,k,im,jm,km);
	int ipjp, jpkp, ipkp; ms.ppx(i,j,k,ipjp,jpkp,ipkp);
	int ipjm, jpkm, ipkm; ms.pmx(i,j,k,ipjm,jpkm,ipkm);
	int imjp, jmkp, imkp; ms.mpx(i,j,k,imjp,jmkp,imkp);
	int imjm, jmkm, imkm; ms.mmx(i,j,k,imjm,jmkm,imkm);

	return 4./9  *   q[id]
	     + 1./9  * ( q[im]   + q[ip]   + q[km]   + q[kp] )
		 + 1./36 * ( q[imkm] + q[ipkp] + q[imkp] + q[ipkm] );

	// return 0.25   *   q[id]
	//      + 0.125  * ( q[im]   + q[ip]   + q[km]   + q[kp] )
	// 	 + 0.0625 * ( q[imkm] + q[ipkp] + q[imkp] + q[ipkm] );
}





// #define DEBUG_VALIDATE_
// #define DEBUG_EFFICIENCY_


#ifdef DEBUG_EFFICIENCY_
#include <sys/time.h>
#include <omp.h>
#include "Bcond.h"


struct timeval *time0 = new struct timeval;
struct timeval *time1 = new struct timeval;
long duration = 0;

int main()
{
	omp_set_num_threads(16);

	int Nx = 33; double Lx = 4;
	int Ny = 33; double Ly = 2;
	int Nz = 33; double Lz = 4;

	Geometry_prdxz geo(Nx,Ny,Nz, Lx,Ly,Lz);
	Mesh ms(geo);
	Scla q(ms), q1(ms);

	geo.InitMesh(0);
	geo.InitInterval();
	geo.InitWaveNumber();
	geo.InitIndices();

	for (int j=0; j<=Ny; j++) {
	for (int k=0; k<=Nz; k++) {
	for (int i=0; i<=Nx; i++) {
		q(i,j,k) = i*j*k;
	}}}

	Bcond::SetBoundaryX(q, 4);
	Bcond::SetBoundaryZ(q, 4);


	gettimeofday(time0, NULL);
	
	for (int j=0; j<=Ny; j++) { cout << j << endl;
	for (int k=0; k<=Nz; k++) {
	for (int i=0; i<=Nx; i++) {
		q1(i,j,k) = Filter::FilterNode(
			ms.xc(i), ms.yc(j), ms.zc(k),
			8*ms.dx(1), 8*ms.dy(1), 8*ms.dz(1),
			q, 0);
	}}}
	
	gettimeofday(time1, NULL);
	duration = 1e6 * (time1->tv_sec - time0->tv_sec) + (time1->tv_usec - time0->tv_usec);

	cout << duration << endl;

	q1.debug_AsciiOutput("", "q1", 0, Ny+1);

	return 0;
}
#endif


#ifdef DEBUG_VALIDATE_
#include "Bcond.h"

int main()
{

	int Nx = 5; double Lx = 4;
	int Ny = 3; double Ly = 2;
	int Nz = 7; double Lz = 4;

	Geometry_prdxz geo(Nx,Ny,Nz, Lx,Ly,Lz);
	// Geometry geo(Nx,Ny,Nz, Lx,Ly,Lz);
	Mesh ms(geo);
	Scla q(ms), q1(ms), q2(ms);

	geo.InitMesh(0);
	geo.InitInterval();
	geo.InitWaveNumber();
	geo.InitIndices();

	// cout << endl << "x, xc" << endl;
	// for (int i=0; i<=Nx; i++) cout << ms.x(i) << ",\t" << ms.xc(i) << endl;
	// cout << endl << "y, yc" << endl;
	// for (int j=0; j<=Ny; j++) cout << ms.y(j) << ",\t" << ms.yc(j) << endl;
	// cout << endl << "z, zc" << endl;
	// for (int k=0; k<=Nz; k++) cout << ms.z(k) << ",\t" << ms.zc(k) << endl;
	// cout << endl;

	for (int j=0; j<=Ny; j++) {
	for (int k=0; k<=Nz; k++) {
	for (int i=0; i<=Nx; i++) {
		q(i,j,k) = i*j*k;
		// q(i,j,k) = i%2 + k%2;
	}}}

	Bcond::SetBoundaryX(q, 4);
	Bcond::SetBoundaryZ(q, 4);

	// Filter::FilterNode(
	// 	ms.xc(1), ms.yc(1), ms.zc(3), 0,0,0, q, 0);

	for (int j=0; j<=Ny; j++) {
	for (int k=0; k<=Nz; k++) {
	for (int i=0; i<=Nx; i++) {
		q1(i,j,k) = Filter::FilterNode(
			ms.xc(i), ms.yc(j), ms.zc(k), 2*ms.dx(1), 0, 2*ms.dz(1),
			// ms.xc(i), ms.yc(j), ms.zc(k), 0,0,0,
			q, 0);

		q2(i,j,k) = Filter::TestFilter(i,j,k,q); // should be 1/4, 2/1, 1/4
	}}}

	q.debug_AsciiOutput("", "q", 0, Ny+1);
	q1.debug_AsciiOutput("", "q1", 0, Ny+1);
	q2.debug_AsciiOutput("", "q2", 0, Ny+1);

	return 0;
}

#endif







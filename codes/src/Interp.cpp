#include "Interp.h"

using namespace std;


double Interp::Shift(double x, double x0, double x1)
// shift x into range [x0,x1)
{
	double l = x1 - x0;
	return x==x1 ? x : fmod(fmod(x-x0, l) + l, l) + x0;
}

double Interp::BiSearch(double x, const double *xs, int i0, int in)
// binary search in range [i0,in) for index of the largest element <= x
{
	int mid, lo = i0, hi = in;
	while (hi-lo > 1)
		if (x < xs[mid=(lo+hi)/2])
			hi = mid;
		else
			lo = mid;
	return lo;
}

double Interp::InterpCell(const Scla &q,
	double x, double xm, double xp, int im, int ip,
	double y, double ym, double yp, int jm, int jp,
	double z, double zm, double zp, int km, int kp )
{
	double dx = xp-xm, dxm = x-xm, dxp = xp-x;
	double dy = yp-ym, dym = y-ym, dyp = yp-y;
	double dz = zp-zm, dzm = z-zm, dzp = zp-z;

	return 1. / (dx * dy * dz) * (
		dxp * dyp * dzp * q(im,jm,km) +
		dxp * dyp * dzm * q(im,jm,kp) +
		dxp * dym * dzp * q(im,jp,km) +
		dxp * dym * dzm * q(im,jp,kp) +
		dxm * dyp * dzp * q(ip,jm,km) +
		dxm * dyp * dzm * q(ip,jm,kp) +
		dxm * dym * dzp * q(ip,jp,km) +
		dxm * dym * dzm * q(ip,jp,kp) );
}

double Interp::InterpNode(double x_, double y_, double z_, const Scla &q, int stgtyp)
{
	const Mesh &ms = q.ms;

	// configure the staggered coordinate to be interpolated
	const double *corx = stgtyp==1 ? ms.x() : ms.xc();
	const double *cory = stgtyp==2 ? ms.y() : ms.yc();
	const double *corz = stgtyp==3 ? ms.z() : ms.zc();

	int i0 = stgtyp==1, in = ms.Nx;
	int j0 = stgtyp==2, jn = ms.Ny;
	int k0 = stgtyp==3, kn = ms.Nz;

	// shift target position into meaningful range
	double x = Shift(x_, ms.x(1), ms.x(ms.Nx));
	double y = Shift(y_, ms.y(1), ms.y(ms.Ny));
	double z = Shift(z_, ms.z(1), ms.z(ms.Nz));
	// double x = fmod(fmod(x_-ms.x(1), ms.Lx) + ms.Lx, ms.Lx) + ms.x(1);
	// double y = fmod(fmod(y_-ms.y(1), ms.Ly) + ms.Ly, ms.Ly) + ms.y(1);
	// double z = fmod(fmod(z_-ms.z(1), ms.Lz) + ms.Lz, ms.Lz) + ms.z(1);

	// binary search for the interpolation range
	int im = BiSearch(x, corx, i0, in), ip = im + 1;
	int jm = BiSearch(y, cory, j0, jn), jp = jm + 1;
	int km = BiSearch(z, corz, k0, kn), kp = km + 1;

	// int mid; // note: im denote the largest point that is <= x
	// while (ip-im > 1) if (x < corx[mid=(im+ip)/2]) ip = mid; else im = mid;
	// while (jp-jm > 1) if (y < cory[mid=(jm+jp)/2]) jp = mid; else jm = mid;
	// while (kp-km > 1) if (z < corz[mid=(km+kp)/2]) kp = mid; else km = mid;

	return InterpCell(q,
		x, corx[im], corx[ip], im, ip,
		y, cory[jm], cory[jp], jm, jp,
		z, corz[km], corz[kp], km, kp );
}

void Interp::InterpBulk(Scla &dst, const Scla &src, int stgtyp)
// whole-bulk linear-interpolation up to virtual boundary
{
	const Mesh &ms = dst.ms;

	// configure the staggered coordinate to be interpolated
	const double *corx = stgtyp==1 ? ms.x() : ms.xc();
	const double *cory = stgtyp==2 ? ms.y() : ms.yc();
	const double *corz = stgtyp==3 ? ms.z() : ms.zc();

	for (int j=(stgtyp==2); j<=ms.Ny; j++) {
	for (int k=(stgtyp==3); k<=ms.Nz; k++) {
	for (int i=(stgtyp==1); i<=ms.Nx; i++) {
		dst(i,j,k) = InterpNode(corx[i], cory[j], corz[k], src, stgtyp);
	}}}
}







// Interp::Interp():
// ms0(ms0),
// ms1(ms1)
// {
// 	isu = new int[ms1.Nx+1]; isa = new int[ms1.Nx+1];
// 	jsv = new int[ms1.Ny+1]; jsa = new int[ms1.Ny+1];
// 	ksw = new int[ms1.Nz+1]; ksa = new int[ms1.Nz+1];
// }

// Interp::~Interp()
// {
// 	delete[] isu; delete[] isa;
// 	delete[] jsv; delete[] jsa;
// 	delete[] ksw; delete[] ksa;
// }


// Interp::InitIndices()
// /* decided interpolation range for ms1 <= ms0 */
// {
// 	int i0, j0, k0;
// 	int i1, j1, k1;

// 	// for every x1, find the largest x0 <= x1
// 	for (i1=i0=1; i1<=ms1.Nx; i1++) {
// 		while (ms0.x(++i0) - ms1.x(i1) < INFTSM)
// 			if (i0 > ms0.Nx) break; // Nx -> Nx
// 		isu[i1] = --i0;
// 	}
// 	for (i1=(i0=0)+1; i1<ms1.Nx; i1++) {
// 		while (ms0.xc(++i0) - ms1.xc(i1) < INFTSM)
// 			if (i0 > ms0.Nx) break; // Nx+1 -> Nx-1
// 		isa[i1] = --i0;
// 	}
// 	// for every y1, find the largest y0 <= y1
// 	for (j1=j0=1; j1<=ms1.Ny; j1++) {
// 		while (ms0.y(++j0) - ms1.y(j1) < INFTSM)
// 			if (j0 > ms0.Ny) break;
// 		jsv[j1] = --j0;
// 	}
// 	for (j1=(j0=0)+1; j1<ms1.Ny; j1++) {
// 		while (ms0.yc(++j0) - ms1.yc(j1) < INFTSM)
// 			if (j0 > ms0.Ny) break;
// 		jsa[j1] = --j0;
// 	}
// 	// for every z1, find the largest z0 <= z1
// 	for (k1=k0=1; k1<=ms1.Nz; k1++) {
// 		while (ms0.z(++k0) - ms1.z(k1) < INFTSM)
// 			if (k0 > ms0.Nz) break;
// 		ksw[k1] = --k0;
// 	}
// 	for (k1=(k0=0)+1; k1<ms1.Nz; k1++) {
// 		while (ms0.zc(++k0) - ms1.zc(k1) < INFTSM)
// 			if (k0 > ms0.Nz) break;
// 		ksa[k1] = --k0;
// 	}

// 	// periodic extrapolation for out-of-range points
// 	// TODO ...
// }





// double Interp::InterpNodeA(double x_, double y_, double z_, const Scla &src)
// {
// 	const Mesh &ms = src.ms;

// 	// shift target position into meaningful range
// 	double x = fmod(fmod(x_, ms.Lx) + ms.Lx, ms.Lx) + ms.x(1);
// 	double y = fmod(fmod(y_, ms.Ly) + ms.Ly, ms.Ly) + ms.y(1);
// 	double z = fmod(fmod(z_, ms.Lz) + ms.Lz, ms.Lz) + ms.z(1);

// 	// binary search for the interpolation range
// 	int im = 0, ip = ms.Nx;
// 	int jm = 0, jp = ms.Ny;
// 	int km = 0, kp = ms.Nz;
// 	int mid;
// 	// note: im denote the largest point that is <= x
// 	while (ip-im > 1) *(x < ms.xc(mid=(im+ip)/2) ? &ip : &im) = mid;
// 	while (jp-jm > 1) *(y < ms.yc(mid=(jm+jp)/2) ? &jp : &jm) = mid;
// 	while (kp-km > 1) *(z < ms.zc(mid=(km+kp)/2) ? &kp : &km) = mid;

// 	return InterpCell(src,
// 		x, ms.xc(im), ms.xc(ip), im, ip,
// 		y, ms.yc(jm), ms.yc(jp), jm, jp,
// 		z, ms.zc(km), ms.zc(kp), km, kp );





// 	// while (ms.xc(++im) - xc < INFTSM) if (im > ms.Nx) break;
// 	// while (ms.yc(++jm) - yc < INFTSM) if (jm > ms.Ny) break;
// 	// while (ms.zc(++km) - zc < INFTSM) if (km > ms.Nz) break;

// 	// double xm = ms.xc(--im), xp = ms.xc(ip=ms.ipa(im));
// 	// double ym = ms.yc(--jm), yp = ms.yc(jp=ms.jpa(jm));
// 	// double zm = ms.zc(--km), zp = ms.zc(kp=ms.kpa(km));

// 	// return InterpCell(src, ms0,
// 	// 	xc,xm,xp,im,ip, yc,ym,yp,jm,jp, zc,zm,zp,km,kp);
// }



// void Interp::Interpolate(Scla &dst, const Scla &src)
// {
// 	int i, im, ip;
// 	int j, jm, jp;
// 	int k, km, kp;
// 	double x, xm, xp, dx, dxm, dxp;
// 	double y, ym, yp, dy, dym, dyp;
// 	double z, zm, zp, dz, dzm, dzp;

// 	for (j=1; j<ms1.Ny; j++) { jm = jsa[j]; jp = ms0.jpa(jm);
// 	for (k=1; k<ms1.Nz; k++) { km = ksa[k]; kp = ms0.kpa(km);
// 	for (i=1; i<ms1.Nx; i++) { im = isa[i]; ip = ms0.ipa(im);

// 		x = ms1.xc(i); xm = ms0.xc(im); xp = ms0.xc(ip);
// 		y = ms1.yc(j); ym = ms0.yc(jm); yp = ms0.yc(jp);
// 		z = ms1.zc(k); zm = ms0.zc(km); zp = ms0.zc(kp);

// 		dx = xp-xm; dxm = x-xm; dxp = xp-x;
// 		dy = yp-ym; dym = y-ym; dyp = yp-y;
// 		dz = zp-zm; dzm = z-zm; dzp = zp-z;

// 		dst(i,j,k) = 1. / (dx * dy * dz) * (
// 			dxp * dyp * dzp * src(im,jm,km) +
// 			dxp * dyp * dzm * src(im,jm,kp) +
// 			dxp * dym * dzp * src(im,jp,km) +
// 			dxp * dym * dzm * src(im,jp,kp) +
// 			dxm * dyp * dzp * src(ip,jm,km) +
// 			dxm * dyp * dzm * src(ip,jm,kp) +
// 			dxm * dym * dzp * src(ip,jp,km) +
// 			dxm * dym * dzm * src(ip,jp,kp) );
// 	}}}
// }



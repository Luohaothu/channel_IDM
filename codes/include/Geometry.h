#pragma once

#include <iostream>
#include <algorithm>
#include <cmath>


class Geometry
{
public:
	const int    Nx, Ny, Nz;
	const double Lx, Ly, Lz;

	const int Nxc, Nzc;
	const int Nxr, Nzr;

	Geometry(int Nx, int Ny, int Nz, double Lx, double Ly, double Lz);
	~Geometry();
	Geometry(const Geometry &geo);
	Geometry& operator=(const Geometry &geo);
	
	void InitMesh(double dy_min, const char *path=NULL);
	void InitInterval();
	void InitIndices();

	void WriteMeshY(const char *path) const;
	void WriteMesh(const char *path) const;

	double *x, *xc;
	double *y, *yc;
	double *z, *zc;

	double *dx, *hx;
	double *dy, *hy;
	double *dz, *hz;

	double *kx, *kx2;
	double *kz, *kz2;

	int *ima, *ipa;
	int *jma, *jpa;
	int *kma, *kpa;

private:
	void DeepCopy(const Geometry &geo);
	void InitMeshY(const char *path);
	void InitMeshY(double dy_min);
	void AlignBoundaryYc(const Geometry &geo);
};


class Geometry_prdz: public Geometry
{
public:
	using Geometry::Geometry;
	void InitMesh(double dy_min, const char *path=NULL);
	void InitIndices();
	void InitWaveNumber();
};


class Geometry_prdxz: public Geometry
{
public:
	using Geometry::Geometry;
	void InitMesh(double dy_min, const char *path=NULL);
	void InitIndices();
	void InitWaveNumber();
};



class Mesh
{
public:
	const int    Nx, Ny, Nz;
	const double Lx, Ly, Lz;

	const int Nxc, Nzc;
	const int Nxr, Nzr;

	Mesh(const Geometry &geo);

	// indexing system
	int idx(int i, int j, int k) const { return j*(Nx+1)*(Nz+1) + k*(Nx+1) + i; };
	int idfxz(int i, int j, int k) const { return idx(0,j,0) + k*Nxr + i; };

	int ima(int i) const { return geo.ima[i]; };
	int jma(int j) const { return geo.jma[j]; };
	int kma(int k) const { return geo.kma[k]; };
	int ipa(int i) const { return geo.ipa[i]; };
	int jpa(int j) const { return geo.jpa[j]; };
	int kpa(int k) const { return geo.kpa[k]; };
	
	// geometry accessing
	double x (int i) const { return geo.x[i]; };
	double y (int j) const { return geo.y[j]; };
	double z (int k) const { return geo.z[k]; };

	double xc(int i) const { return geo.xc[i]; };
	double yc(int j) const { return geo.yc[j]; };
	double zc(int k) const { return geo.zc[k]; };

	const double* x() const { return geo.x; };
	const double* y() const { return geo.y; };
	const double* z() const { return geo.z; };
	const double* xc() const { return geo.xc; };
	const double* yc() const { return geo.yc; };
	const double* zc() const { return geo.zc; };

	double hx(int i) const { return geo.hx[i]; };
	double hy(int j) const { return geo.hy[j]; };
	double hz(int k) const { return geo.hz[k]; };
	double dx(int i) const { return geo.dx[i]; };
	double dy(int j) const { return geo.dy[j]; };
	double dz(int k) const { return geo.dz[k]; };

	double hx(int i, double &hxm, double &hxp) const { hxm = hx(i-1); hxp = hx(i+1); return hx(i); };
	double hy(int j, double &hym, double &hyp) const { hym = hy(j-1); hyp = hy(j+1); return hy(j); };
	double hz(int k, double &hzm, double &hzp) const { hzm = hz(k-1); hzp = hz(k+1); return hz(k); };
	double dx(int i, double &dxm, double &dxp) const { dxm = dx(i-1); dxp = dx(i+1); return dx(i); };
	double dy(int j, double &dym, double &dyp) const { dym = dy(j-1); dyp = dy(j+1); return dy(j); };
	double dz(int k, double &dzm, double &dzp) const { dzm = dz(k-1); dzp = dz(k+1); return dz(k); };
	
	double kx (int i) const { return geo.kx [i]; };
	double kz (int k) const { return geo.kz [k]; };
	double kx2(int i) const { return geo.kx2[i]; };
	double kz2(int k) const { return geo.kz2[k]; };

	// batch indexing
	void ipx(int i, int j, int k, int &ip, int &jp, int &kp) const;
	void imx(int i, int j, int k, int &im, int &jm, int &km) const;

	void ppx(int i, int j, int k, int &ipjp, int &jpkp, int &ipkp) const;
	void pmx(int i, int j, int k, int &ipjm, int &jpkm, int &ipkm) const;
	void mpx(int i, int j, int k, int &imjp, int &jmkp, int &imkp) const;
	void mmx(int i, int j, int k, int &imjm, int &jmkm, int &imkm) const;

	void dcx(int i, int j, int k, double &dxc, double &dyc, double &dzc) const;
	void dpx(int i, int j, int k, double &dxp, double &dyp, double &dzp) const;
	void dmx(int i, int j, int k, double &dxm, double &dym, double &dzm) const;

	void hcx(int i, int j, int k, double &hxc, double &hyc, double &hzc) const;
	void hpx(int i, int j, int k, double &hxp, double &hyp, double &hzp) const;
	void hmx(int i, int j, int k, double &hxm, double &hym, double &hzm) const;

private:
	const Geometry geo;

};





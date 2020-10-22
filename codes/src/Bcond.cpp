#include "Bcond.h"
#include "Interp.h"
#include "Filter.h"

using namespace std;


double SlipLength(int i, int j, int k, const Vctr &vel, const Flow &vis, const Vctr &shear);
void ShearStress(Vctr &shear, const Vctr &vel, const Flow &vis);


// ***** configure boundary conditions ***** //

void Bcond::ChannelNoSlip(Boundaries &bc, Boundaries &sbc, const Mesh &ms)
{
	double dy0 = ms.dy(0), dyn = ms.dy(ms.Ny);
	double dy1 = ms.dy(1), dym = ms.dy(ms.Ny-1);

	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		bc.ub3(i,k) = bc.vb3(i,k) = bc.wb3(i,k) = 0;
		bc.ub4(i,k) = bc.vb4(i,k) = bc.wb4(i,k) = 0;

		sbc.ub3(i,k) = sbc.wb3(i,k) = dy0 / dy1;
		sbc.ub4(i,k) = sbc.wb4(i,k) = dyn / dym;

		sbc.vb3(i,k) = 0;
		sbc.vb4(i,k) = 0;
	}}
}

void Bcond::ChannelDirichlet(Boundaries &bc, Boundaries &sbc, const Mesh &ms, const Vctr &vel)
{
	const Mesh &ms0 = vel.ms;

	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	#pragma omp parallel for
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		sbc.ub3(i,k) = sbc.vb3(i,k) = sbc.wb3(i,k) = 0;
		sbc.ub4(i,k) = sbc.vb4(i,k) = sbc.wb4(i,k) = 0;

		// // directly assign off-wall boundary. note: must be staggered
		// int jb1 =(ms0.Ny - ms.Ny) / 2;
		// int jb2 = ms0.Ny - jb1;

		// bc.ub3(i,k) = u(i,jb1,  k);
		// bc.vb3(i,k) = v(i,jb1+1,k);
		// bc.wb3(i,k) = w(i,jb1,  k);

		// bc.ub4(i,k) = u(i,jb2,k);
		// bc.vb4(i,k) = v(i,jb2,k);
		// bc.wb4(i,k) = w(i,jb2,k);

		// interpolate off-wall boundary
		bc.ub3(i,k) = Interp::InterpNodeU(ms.x(i),ms.yc(0),ms.zc(k), u);
		bc.vb3(i,k) = Interp::InterpNodeV(ms.xc(i),ms.y(1),ms.zc(k), v);
		bc.wb3(i,k) = Interp::InterpNodeW(ms.xc(i),ms.yc(0),ms.z(k), w);

		bc.ub4(i,k) = Interp::InterpNodeU(ms.x(i),ms.yc(ms.Ny),ms.zc(k), u);
		bc.vb4(i,k) = Interp::InterpNodeV(ms.xc(i),ms.y(ms.Ny),ms.zc(k), v);
		bc.wb4(i,k) = Interp::InterpNodeW(ms.xc(i),ms.yc(ms.Ny),ms.z(k), w);
	}}
}

void Bcond::ChannelDirichlet(Boundaries &bc, Boundaries &sbc, const Mesh &ms, const Vctr &vel, double rsclx, double rsclu)
{
	const Mesh &ms0 = vel.ms;

	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	#pragma omp parallel for
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		sbc.ub3(i,k) = sbc.vb3(i,k) = sbc.wb3(i,k) = 0;
		sbc.ub4(i,k) = sbc.vb4(i,k) = sbc.wb4(i,k) = 0;

		double x = ms.x(i) * rsclx, xc = ms.xc(i) * rsclx;
		double z = ms.z(k) * rsclx, zc = ms.zc(k) * rsclx;
		double dx= ms.dx(i)* rsclx, hx = ms.hx(i) * rsclx;
		double dz= ms.dz(k)* rsclx, hz = ms.hz(k) * rsclx;

		double y  = Filter::WallRscl(ms.y (1), rsclx);
		double yc = Filter::WallRscl(ms.yc(0), rsclx);

		// filter off-wall boundary
		bc.ub3(i,k) = Filter::FilterNodeU(x,yc,zc, hx,0,dz, u) * rsclu;
		bc.vb3(i,k) = Filter::FilterNodeV(xc,y,zc, dx,0,dz, v) * rsclu;
		bc.wb3(i,k) = Filter::FilterNodeW(xc,yc,z, dx,0,hz, w) * rsclu;

		y  = Filter::WallRscl(ms.y (ms.Ny), rsclx);
		yc = Filter::WallRscl(ms.yc(ms.Ny), rsclx);

		bc.ub4(i,k) = Filter::FilterNodeU(x,yc,zc, hx,0,dz, u) * rsclu;
		bc.vb4(i,k) = Filter::FilterNodeV(xc,y,zc, dx,0,dz, v) * rsclu;
		bc.wb4(i,k) = Filter::FilterNodeW(xc,yc,z, dx,0,hz, w) * rsclu;

		// handle HALFMFU
		if (ms0.Ly < 2) bc.vb4(i,k) *= -1;
	}}
}

void Bcond::ChannelRobin(Boundaries &bc, Boundaries &sbc, const Mesh &ms)
{
	double dy0 = ms.dy(0), dyn = ms.dy(ms.Ny);
	double dy1 = ms.dy(1), dym = ms.dy(ms.Ny-1);

	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		double lu = 0.008;
		double lv = 0.008;
		double lw = 0.008;

		bc.ub3(i,k) = bc.vb3(i,k) = bc.wb3(i,k) = 0;
		bc.ub4(i,k) = bc.vb4(i,k) = bc.wb4(i,k) = 0;

		sbc.ub3(i,k) = sbc.wb3(i,k) = (dy0 - 2*lu) / (dy1 + 2*lu);
		sbc.ub4(i,k) = sbc.wb4(i,k) = (dyn - 2*lu) / (dym + 2*lu);

		sbc.vb3(i,k) = - lv / (dy1 + lv);
		sbc.vb4(i,k) = - lv / (dym + lv);
	}}
}

void Bcond::ChannelRobin(Boundaries &bc, Boundaries &sbc,
	const Vctr &vel, const Flow &vis, const Vctr &vel0, const Flow &vis0)
{
	const Mesh &ms = vel.ms;

	// solve for reference shear stress from DNS data
	Vctr shear0(vel0.ms);
	ShearStress(shear0, vel0, vis0);

	for (int j=0; j<=ms.Ny; j++) {
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {
		shear0[1](i,j,k) = j < ms.Ny/2 ? 0.004 : -0.004;
		shear0[2](i,j,k) = 0;
		shear0[3](i,j,k) = 0;
	}}} // TODO:  test this


	// solve for slip length and corresponding BC
	double dy0 = ms.dy(0), dyn = ms.dy(ms.Ny);
	double dy1 = ms.dy(1), dym = ms.dy(ms.Ny-1);

	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		bc.ub3(i,k) = bc.vb3(i,k) = bc.wb3(i,k) = 0;
		bc.ub4(i,k) = bc.vb4(i,k) = bc.wb4(i,k) = 0;

		double l_bot = SlipLength(i,1,    k, vel, vis, shear0);//0.008;//
		double l_top = SlipLength(i,ms.Ny,k, vel, vis, shear0);//0.008;//

		if (i==10 && k == 10) cout << l_bot << '\t' << l_top << endl;

		sbc.ub3(i,k) = sbc.wb3(i,k) = (dy0 - 2*l_bot) / (dy1 + 2*l_bot);
		sbc.ub4(i,k) = sbc.wb4(i,k) = (dyn - 2*l_top) / (dym + 2*l_top);

		sbc.vb3(i,k) = - l_bot / (dy1 + l_bot);
		sbc.vb4(i,k) = - l_top / (dym + l_top);
	}}
}

void Bcond::ChannelHalf(Boundaries &bc, Boundaries &sbc, const Mesh &ms)
{
	double dy0 = ms.dy(0);
	double dy1 = ms.dy(1);

	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		bc.ub3(i,k) = 0; sbc.ub3(i,k) = dy0 / dy1; // for points deployed on virtual boundary,
		bc.wb3(i,k) = 0; sbc.wb3(i,k) = dy0 / dy1; // interpolation is needed to satisfy BC on real boundary
		bc.vb3(i,k) = 0; sbc.vb3(i,k) = 0;

		bc.ub4(i,k) = 0; sbc.ub4(i,k) = -1;
		bc.wb4(i,k) = 0; sbc.wb4(i,k) = -1;
		bc.vb4(i,k) = 0; sbc.vb4(i,k) = 0;
	}}
}

void Bcond::TblCycling(Boundaries &bc, Boundaries &sbc, const Mesh &ms, double ufree)
{
	double dy0 = ms.dy(0);
	double dy1 = ms.dy(1);

	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		bc.ub3(i,k) = 0; sbc.ub3(i,k) = dy0 / dy1;
		bc.wb3(i,k) = 0; sbc.wb3(i,k) = dy0 / dy1;
		bc.vb3(i,k) = 0; sbc.vb3(i,k) = 0;

		bc.ub4(i,k) = ufree; sbc.ub4(i,k) = 0;
		bc.vb4(i,k) = 0; sbc.vb4(i,k) = -1;
		bc.wb4(i,k) = 0; sbc.wb4(i,k) = -1;
	}}
}


// ***** apply boundary conditions for velocity fields ***** //

void Bcond::SetBoundaryX(Vctr &vel, const Boundaries &bc, const Boundaries &sbc)
{
	const Mesh &ms = vel.ms;

	int Nx = ms.Nx;
	int Ny = ms.Ny;
	int Nz = ms.Nz;

	Scla &u = vel[1];
	Scla &v = vel[2];
	Scla &w = vel[3];

	for (int j=0; j<=Ny; j++) {
	for (int k=0; k<=Nz; k++) {

		double ubc1 = bc.ub1(j,k), sbu1 = sbc.ub1(j,k);
		double vbc1 = bc.vb1(j,k), sbv1 = sbc.vb1(j,k);
		double wbc1 = bc.wb1(j,k), sbw1 = sbc.wb1(j,k);

		double ubc2 = bc.ub2(j,k), sbu2 = sbc.ub2(j,k);
		double vbc2 = bc.vb2(j,k), sbv2 = sbc.vb2(j,k);
		double wbc2 = bc.wb2(j,k), sbw2 = sbc.wb2(j,k);

		u(0,j,k) = 0;

		u(1,j,k) = ubc1 - sbu1 * u(2,j,k);
		v(0,j,k) = vbc1 - sbv1 * v(1,j,k);
		w(0,j,k) = wbc1 - sbw1 * w(1,j,k);

		u(Nx,j,k) = ubc2 - sbu2 * u(Nx-1,j,k);
		v(Nx,j,k) = vbc2 - sbv2 * v(Nx-1,j,k);
		w(Nx,j,k) = wbc2 - sbw2 * w(Nx-1,j,k);
	}}
}

void Bcond::SetBoundaryY(Vctr &vel, const Boundaries &bc, const Boundaries &sbc)
{
	const Mesh &ms = vel.ms;

	int Nx = ms.Nx;
	int Ny = ms.Ny;
	int Nz = ms.Nz;

	Scla &u = vel[1];
	Scla &v = vel[2];
	Scla &w = vel[3];

	for (int k=0; k<=Nz; k++) {
	for (int i=0; i<=Nx; i++) {

		double ubc3 = bc.ub3(i,k), sbu3 = sbc.ub3(i,k);
		double vbc3 = bc.vb3(i,k), sbv3 = sbc.vb3(i,k);
		double wbc3 = bc.wb3(i,k), sbw3 = sbc.wb3(i,k);

		double ubc4 = bc.ub4(i,k), sbu4 = sbc.ub4(i,k);
		double vbc4 = bc.vb4(i,k), sbv4 = sbc.vb4(i,k);
		double wbc4 = bc.wb4(i,k), sbw4 = sbc.wb4(i,k);

		v(i,0,k) = 0;

		u(i,0,k) = ubc3 - sbu3 * u(i,1,k);
		v(i,1,k) = vbc3 - sbv3 * v(i,2,k);
		w(i,0,k) = wbc3 - sbw3 * w(i,1,k);

		u(i,Ny,k) = ubc4 - sbu4 * u(i,Ny-1,k);
		v(i,Ny,k) = vbc4 - sbv4 * v(i,Ny-1,k);
		w(i,Ny,k) = wbc4 - sbw4 * w(i,Ny-1,k);
	}}
}


void Bcond::SetBoundaryX(Vctr &vel)
{
	const Mesh &ms = vel.ms;

	int Nx = ms.Nx;
	int Ny = ms.Ny;
	int Nz = ms.Nz;

	Scla &u = vel[1];
	Scla &v = vel[2];
	Scla &w = vel[3];

	for (int k=0; k<=Nz; k++) {
	for (int j=0; j<=Ny; j++) {

		u(0,j,k) = 0;
		v(0,j,k) = v(Nx-1,j,k);
		w(0,j,k) = w(Nx-1,j,k);

		u(Nx,j,k) = u(1,j,k);
		v(Nx,j,k) = v(1,j,k);
		w(Nx,j,k) = w(1,j,k);
	}}
}

void Bcond::SetBoundaryZ(Vctr &vel)
{
	const Mesh &ms = vel.ms;

	int Nx = ms.Nx;
	int Ny = ms.Ny;
	int Nz = ms.Nz;

	Scla &u = vel[1];
	Scla &v = vel[2];
	Scla &w = vel[3];

	for (int j=0; j<=Ny; j++) {
	for (int i=0; i<=Nx; i++) {

		u(i,j,0) = u(i,j,Nz-1);
		v(i,j,0) = v(i,j,Nz-1);
		w(i,j,0) = 0;

		u(i,j,Nz) = u(i,j,1);
		v(i,j,Nz) = v(i,j,1);
		w(i,j,Nz) = w(i,j,1);
	}}
}


// ***** set homogeneous boundary for cell-centered scalars ***** //

void Bcond::SetBoundaryX(Scla &q, int ord)
{
	const Mesh &ms = q.ms;

	double h1 = ms.hx(1);
	double h2 = ms.hx(2);
	double hm = ms.hx(ms.Nx-1);
	double hn = ms.hx(ms.Nx);

	for (int j=0; j<=ms.Ny; j++) {
	for (int k=0; k<=ms.Nz; k++) {

		int i0 = ms.idx(0,j,k), in = ms.idx(ms.Nx,j,k);
		int i1 = ms.idx(1,j,k), im = ms.idx(ms.Nx-1,j,k);
		int i2 = ms.idx(2,j,k), il = ms.idx(ms.Nx-2,j,k);

		// homogeneous ord-th order derivative or periodic
		q[i0] = ord==0 ? 0 : ord==1 ? q[i1] : ord==2 ? (h1/h2+1) * q[i1] - h1/h2 * q[i2] : q[im];
		q[in] = ord==0 ? 0 : ord==1 ? q[im] : ord==2 ? (hn/hm+1) * q[im] - hn/hm * q[il] : q[i1];
	}}
}

void Bcond::SetBoundaryY(Scla &q, int ord)
{
	const Mesh &ms = q.ms;

	double h1 = ms.hy(1);
	double h2 = ms.hy(2);
	double hm = ms.hy(ms.Ny-1);
	double hn = ms.hy(ms.Ny);

	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		int j0 = ms.idx(i,0,k), jn = ms.idx(i,ms.Ny,k);
		int j1 = ms.idx(i,1,k), jm = ms.idx(i,ms.Ny-1,k);
		int j2 = ms.idx(i,2,k), jl = ms.idx(i,ms.Ny-2,k);

		// homogeneous ord-th order derivative or periodic
		q[j0] = ord==0 ? 0 : ord==1 ? q[j1] : ord==2 ? (h1/h2+1) * q[j1] - h1/h2 * q[j2] : q[jm];
		q[jn] = ord==0 ? 0 : ord==1 ? q[jm] : ord==2 ? (hn/hm+1) * q[jm] - hn/hm * q[jl] : q[j1];
	}}
}

void Bcond::SetBoundaryZ(Scla &q, int ord)
{
	const Mesh &ms = q.ms;

	double h1 = ms.hz(1);
	double h2 = ms.hz(2);
	double hm = ms.hz(ms.Nz-1);
	double hn = ms.hz(ms.Nz);

	for (int j=0; j<=ms.Ny; j++) {
	for (int i=0; i<=ms.Nx; i++) {

		int k0 = ms.idx(i,j,0), kn = ms.idx(i,j,ms.Nz);
		int k1 = ms.idx(i,j,1), km = ms.idx(i,j,ms.Nz-1);
		int k2 = ms.idx(i,j,2), kl = ms.idx(i,j,ms.Nz-2);

		// homogeneous ord-th order derivative or periodic
		q[k0] = ord==0 ? 0 : ord==1 ? q[k1] : ord==2 ? (h1/h2+1) * q[k1] - h1/h2 * q[k2] : q[km];
		q[kn] = ord==0 ? 0 : ord==1 ? q[km] : ord==2 ? (hn/hm+1) * q[km] - hn/hm * q[kl] : q[k1];
	}}
}








void DeltaTau(double dtau[3], int j, const Vctr &vel, const Flow &vis, const Vctr &shear)
{
	const Mesh &ms = vel.ms; // note: shear has different mesh

	const Scla &u = vel[1], &nux = vis.SeeVec(1), &tau12 = shear[1];
	const Scla &v = vel[2], &nuy = vis.SeeVec(2), &tau23 = shear[2];
	const Scla &w = vel[3], &nuz = vis.SeeVec(3), &tau13 = shear[3];

	dtau[0] = 0;
	dtau[1] = 0;
	dtau[2] = 0;

	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double x = ms.x(i), xc = ms.xc(i), dx = ms.dx(i), hx = ms.hx(i);
		double z = ms.z(k), zc = ms.zc(k), dz = ms.dz(k), hz = ms.hz(k);
		double y = ms.y(j);
		double offsetx = .5 * tau12.ms.dx(1); // filter on edge requires offset
		double offsetz = .5 * tau23.ms.dz(1);

		// reference shear stress filtered from DNS data
		double tau12dns = j==1 ? 0.004 : -0.004;//Filter::FilterNodeV(x+offsetx, y, zc, hx,0,dz, tau12);
		double tau23dns = 0;//Filter::FilterNodeV(xc, y, z+offsetz, dx,0,hz, tau23);

		// strain rate of the current velocity field
		const double *sr = vel.ShearStrain(i,j,k);

		// Fij = <tau^DNS_ij - tau_ij>
		dtau[0] += (tau12dns - 2 * nuz(i,j,k) * sr[0]) / ((ms.Nx-1) * (ms.Nz-1));
		dtau[1] += (tau23dns - 2 * nux(i,j,k) * sr[1]) / ((ms.Nx-1) * (ms.Nz-1));
		// lack Reynolds stress !!!
	}}
}


double SlipLength(int i, int j, int k, const Vctr &vel, const Flow &vis, const Vctr &shear)
// calculate slip length of a cell up to ral boundaries
// based on least-square of stress-like tensor deployed on cell edges
{
	const Mesh &ms = vel.ms; // note: shear has different mesh

	const Scla &u = vel[1], &nux = vis.SeeVec(1), &tau12 = shear[1];
	const Scla &v = vel[2], &nuy = vis.SeeVec(2), &tau23 = shear[2];
	const Scla &w = vel[3], &nuz = vis.SeeVec(3), &tau13 = shear[3];

	int id =              ms.idx(i,j,k);
	int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
	int im, jm, km;       ms.imx(i,j,k,im,jm,km);
	double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
	double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
	double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);

	// derivative on edges
	double dudy = (u[id] - u[jm]) / hyc;
	double dvdx = (v[id] - v[im]) / hxc;
	double dwdy = (w[id] - w[jm]) / hyc;
	double dvdz = (v[id] - v[km]) / hzc;

	// reference shear stress filtered from DNS data
	double x = ms.x(i), xc = ms.xc(i), dx = ms.dx(i), hx = ms.hx(i);
	double z = ms.z(k), zc = ms.zc(k), dz = ms.dz(k), hz = ms.hz(k);
	double y = ms.y(j);
	double offsetx = .5 * tau12.ms.dx(1); // filter on edge requires offset
	double offsetz = .5 * tau23.ms.dz(1);

	double tau12dns = Filter::FilterNodeV(x+offsetx, y, zc, hx,0,dz, tau12);
	double tau23dns = Filter::FilterNodeV(xc, y, z+offsetz, dx,0,hz, tau23);

	// Fij = tau^DNS_ij - tau_ij
	double f12 = tau12dns - nuz[id] * (dvdx + dudy);
	double f23 = tau23dns - nux[id] * (dvdz + dwdy);

	// Lij = u_i * u_j
	double l12 = .25/hxc/hyc * (u[id]*dym + u[jm]*dyc) * (v[id]*dxm + v[im]*dxc);
	double l23 = .25/hzc/hyc * (w[id]*dym + w[jm]*dyc) * (v[id]*dzm + v[km]*dzc);

	// Mij = du_i/dy * du_j/dy
	double m12, m23;
	if (j == 1) {
		m12 = dudy * .5/hxc/dyc * ((v[jp]-v[id]) * dxm + (v(ms.ima(i),ms.jpa(j),k)-v[im]) * dxc);
		m23 = dwdy * .5/hzc/dyc * ((v[jp]-v[id]) * dzm + (v(i,ms.jpa(j),ms.kma(k))-v[km]) * dzc);
	}
	else if (j == ms.Ny) {
		m12 = dudy * .5/hxc/dym * ((v[id]-v[jm]) * dxm + (v[im]-v(ms.ima(i),ms.jma(j),k)) * dxc);
		m23 = dwdy * .5/hzc/dym * ((v[id]-v[jm]) * dzm + (v[km]-v(i,ms.jma(j),ms.kma(k))) * dzc);
	}
	else {
		m12 = dudy * .25/hxc/hyc * ((v[jp]-v[jm]) * dxm + (v(ms.ima(i),ms.jpa(j),k)-v(ms.ima(i),ms.jma(j),k)) * dxc);
		m23 = dwdy * .25/hzc/hyc * ((v[jp]-v[jm]) * dzm + (v(i,ms.jpa(j),ms.kma(k))-v(i,ms.jma(j),ms.kma(k))) * dzc);
	}

	double l_sq = ((l12+f12)*m12 + (l23+f23)*m23) / (m12*m12 + m23*m23);

	// return sqrt(fmin(fmax(l_sq, 0), 1));
	return sqrt(fmin(fabs(l_sq), 1));	
}


#define DIRTY_TRICK_BCOND_

void ShearStress(Vctr &shear, const Vctr &vel, const Flow &vis)
// solve shear stress on edges up to virtual boundaries
{
	const Mesh &ms = shear.ms;
	const Scla &nux = vis.SeeVec(1);
	const Scla &nuy = vis.SeeVec(2);
	const Scla &nuz = vis.SeeVec(3);

	Scla &tau12 = shear[1];
	Scla &tau23 = shear[2];
	Scla &tau13 = shear[3];

#ifndef DIRTY_TRICK_BCOND_
	for (int j=0; j<=ms.Ny; j++) {
#else
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
#endif
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		const double *sr = vel.ShearStrain(i,j,k);

		if (i>0 && j>0) tau12(i,j,k) = 2 * nuz(i,j,k) * sr[0];
		if (j>0 && k>0) tau23(i,j,k) = 2 * nux(i,j,k) * sr[1];
		if (i>0 && k>0) tau13(i,j,k) = 2 * nuy(i,j,k) * sr[2];
	}}}
}



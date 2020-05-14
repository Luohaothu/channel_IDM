#include "Bcond.h"
#include "Interp.h"
#include "Filter.h"


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

		// // filter off-wall boundary
		// bc.ub3(i,k) = Filter::FilterNodeU(ms.x(i),ms.yc(0),ms.zc(k), ms.hx(i),0,ms.dz(k), u);
		// bc.vb3(i,k) = Filter::FilterNodeV(ms.xc(i),ms.y(1),ms.zc(k), ms.dx(i),0,ms.dz(k), v);
		// bc.wb3(i,k) = Filter::FilterNodeW(ms.xc(i),ms.yc(0),ms.z(k), ms.dx(i),0,ms.hz(k), w);

		// bc.ub4(i,k) = Filter::FilterNodeU(ms.x(i),ms.yc(ms.Ny),ms.zc(k), ms.hx(i),0,ms.dz(k), u);
		// bc.vb4(i,k) = Filter::FilterNodeV(ms.xc(i),ms.y(ms.Ny),ms.zc(k), ms.dx(i),0,ms.dz(k), v);
		// bc.wb4(i,k) = Filter::FilterNodeW(ms.xc(i),ms.yc(ms.Ny),ms.z(k), ms.dx(i),0,ms.hz(k), w);
	}}
}

void Bcond::ChannelRobin(Boundaries &bc, Boundaries &sbc, const Mesh &ms, const Vctr &vel)
{
	const Mesh &ms0 = vel.ms;
	
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







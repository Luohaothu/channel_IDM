#include "WM.h"
#include "SGS.h"
#include "Interp.h"
#include "Filter.h"
#include "Bcond.h"

using namespace std;


void WM::UniformWallShear(Flow &vis, const Vctr &vel, double tau12)
// rescale viscosity on z-edge to ensure uniform tau12 on every boundary point
{
	const Mesh &ms = vis.ms;

	Scla &nuz = vis.GetVec(3);

	for (int j=1; j<=ms.Ny; j+= ms.Ny-1) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double s12 = vel.ShearStrain(i,j,k)[0];
		nuz(i,j,k) = fmax(.5*tau12/s12, 0);
	}}}
}


double ReynoldsStressDefect(int j, const Vctr &vel, const Vctr &veldns, double rsclx, double rsclu)
// calculate Reynolds stress R12 - R12^r at Z-edge of layer j
{
	const Mesh &ms= vel.ms, &ms0= veldns.ms;
	const Scla &u = vel[1], &u0 = veldns[1];
	const Scla &v = vel[2], &v0 = veldns[2];

	// the wall normal position on which stresses act
	const double y = ms.y(j);
	const double y0 = WallRscl(y, rsclx); // position in DNS field with y^+ matched

	double r12dns = 0; // DNS Reynolds stress
	double r12les = 0; // LES resolved Reynolds

	// calculate reference (DNS) Reynolds shear stress
	int j3 = Interp::BiSearch(y0, ms0.yc(), 0, ms0.Ny);
	int j4 = j3 + 1;

	double um3 = u0.meanxz(j3), y3 = ms0.yc(j3);
	double um4 = u0.meanxz(j4), y4 = ms0.yc(j4);
	double vm3 = .5 * (v0.meanxz(j3) + v0.meanxz(ms0.jpa(j3)));
	double vm4 = .5 * (v0.meanxz(j4) + v0.meanxz(ms0.jpa(j4))); // note: if j4==ms0.Ny, jpa is not valid

	#pragma omp parallel for reduction(+: r12dns)
	for (int k=1; k<ms0.Nz; k++) { double z = ms0.zc(k);
	for (int i=1; i<ms0.Nx; i++) { double x = ms0.xc(i);
		// calculate uv on cell-centers and interpolate to the desired y position
		// note: interpolation must be afterwards, otherwise the Reynolds stress would be defected
		r12dns += 1. / (y4-y3) * (
			(Interp::InterpNodeU(x,y3,z,u0) - um3)
		  * (Interp::InterpNodeV(x,y3,z,v0) - vm3) * (y4-y0)
		  + (Interp::InterpNodeU(x,y4,z,u0) - um4)
		  * (Interp::InterpNodeV(x,y4,z,v0) - vm4) * (y0-y3) );
	}}

	r12dns *= rsclu * rsclu / (ms0.Nx-1.) / (ms0.Nz-1.);

	// calculate resolved Reynolds stress
	um3 = u.meanxz(j3 = ms.jma(j)); y3 = ms.yc(j3);
	um4 = u.meanxz(j4 = j);         y4 = ms.yc(j4);
	double um = (um3 * (y4-y) + um4 * (y-y3)) / (y4-y3);
	double vm = v.meanxz(j);

	#pragma omp parallel for reduction(+: r12les)
	for (int k=1; k<ms.Nz; k++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hyc=ms.hy(j);
	for (int i=1; i<ms.Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxc=ms.hx(i);

		int id =        ms.idx(i,j,k);
		int im, jm, km; ms.imx(i,j,k,im,jm,km);

		r12les += (.5/hyc * (u[id]*dym + u[jm]*dyc) - um)
		        * (.5/hxc * (v[id]*dxm + v[im]*dxc) - vm);
	}}

	r12les /= (ms.Nx-1.) * (ms.Nz-1.);

	return r12dns - r12les;
}


void ProcessBoundaryVis(Flow &vis)
{
	Scla &nuc = vis.GetScl();
	Bcond::SetBoundaryY(nuc, 1); // homogeneous Neumann
	Bcond::SetBoundaryX(nuc, 3); // periodic
	Bcond::SetBoundaryZ(nuc, 3); // periodic
	vis.CellCenter2Edge();
}


void OffWallSubGridShear(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double rsclx, double rsclu)
// supply sgs stress 12 & 23 on off-wall boundary through viscosity
{
	const Mesh &ms = vis.ms, &ms0 = veldns.ms;

	const Scla &u = vel[1], &u0 = veldns[1];
	const Scla &v = vel[2], &v0 = veldns[2];
	const Scla &w = vel[3], &w0 = veldns[3];

	// viscosity & sgs-stress on edges aimed for
	Vctr shearsgs(ms);
	const Scla &tau23sgs = shearsgs[2];
	const Scla &tau12sgs = shearsgs[1];

	Scla &nux = vis.GetVec(1);
	Scla &nuz = vis.GetVec(3);

	// sgs shear stress filtered from resolved velocity field
	SGS::SubGridShearStress(shearsgs, veldns, rsclx, rsclu);

	// for top & bottom real boundaries
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

		// the wall normal position on which stresses act
		const double y = ms.y(j), y0 = WallRscl(y, rsclx);

		// ***** rescale the DNS-filtered sgs stress by Reynolds stress defect ***** //

		double r12dns = 0; // DNS Reynolds stress
		double r12les = 0; // LES resolved Reynolds
		double s12sgs = 0; // mean sgs-stress filtered from DNS

		// calculate reference (DNS) Reynolds shear stress
		int j3 = Interp::BiSearch(y0, ms0.yc(), 0, ms0.Ny);
		int j4 = j3 + 1;

		double um3 = u0.meanxz(j3), y3 = ms0.yc(j3);
		double um4 = u0.meanxz(j4), y4 = ms0.yc(j4);
		double vm3 = .5 * (v0.meanxz(j3) + v0.meanxz(ms0.jpa(j3)));
		double vm4 = .5 * (v0.meanxz(j4) + v0.meanxz(ms0.jpa(j4))); // note: if j4==ms0.Ny, jpa is not valid

		#pragma omp parallel for reduction(+: r12dns)
		for (int k=1; k<ms0.Nz; k++) { double z = ms0.zc(k);
		for (int i=1; i<ms0.Nx; i++) { double x = ms0.xc(i);
			// calculate uv on cell-centers and interpolate to the desired y position
			// note: interpolation must be afterwards, otherwise the Reynolds stress would be defected
			r12dns += 1. / (y4-y3) * (
				(Interp::InterpNodeU(x,y3,z,u0) - um3)
			  * (Interp::InterpNodeV(x,y3,z,v0) - vm3) * (y4-y0)
			  + (Interp::InterpNodeU(x,y4,z,u0) - um4)
			  * (Interp::InterpNodeV(x,y4,z,v0) - vm4) * (y0-y3) );
		}}

		r12dns *= rsclu * rsclu
			   * (ms.Nx-1.) / (ms0.Nx-1.) // rescale to match the other two
			   * (ms.Nz-1.) / (ms0.Nz-1.);

		// calculate resolved Reynolds stress and mean sgs-stress
		um3 = u.meanxz(j3 = ms.jma(j)); y3 = ms.yc(j3);
		um4 = u.meanxz(j4 = j);         y4 = ms.yc(j4);
		double vm  = v.meanxz(j);
		double um  = (um3 * (y4-y) + um4 * (y-y3)) / (y4-y3);

		for (int k=1; k<ms.Nz; k++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hyc=ms.hy(j);
		for (int i=1; i<ms.Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxc=ms.hx(i);

			int id =        ms.idx(i,j,k);
			int im, jm, km; ms.imx(i,j,k,im,jm,km);

			s12sgs += tau12sgs[id];
			r12les += (.5/hyc * (u[id]*dym + u[jm]*dyc) - um)
			        * (.5/hxc * (v[id]*dxm + v[im]*dxc) - vm);
		}}

		// ***** rescale viscosity to account for low-order differencing error ***** //

		// // Option 1: construct low- & high-order differenced U gradient
		// int j5u = j==1 ? j4u+1 : j3u-1;
		// double um5 = u.meanxz(j5u);
		// double y5u = ms.yc(j5u);

		// double dyU2 = (y4u+y5u-2*y) / (y4u-y3u) / (y3u-y5u) * um3
		//             + (y3u+y5u-2*y) / (y3u-y4u) / (y4u-y5u) * um4
		//             + (y3u+y4u-2*y) / (y3u-y5u) / (y5u-y4u) * um5;

		// // Option 2: rescale physical viscosity by replacing low-order differencing with log law
		// double dyU2 = (um4-um3) / log((1-fabs(y4u-1))/(1-fabs(y3u-1))) / (1-fabs(y-1));
		// double dyU1 = (um4-um3) / (y4u-y3u);
		// const double rsclvis = fabs(dyU2 / dyU1);

		const double rsclsgs = fmin(fabs((r12dns - r12les) / s12sgs), 2.);
		const double rsclvis = fabs((y4-y3) / (1-fabs(y-1)) / log((1-fabs(y4-1))/(1-fabs(y3-1))));


		FILE *fp = fopen("WMLOG.dat", "a");
		if (j==1) fprintf(fp, "%f\t%f\t", rsclsgs, rsclvis);
		else      fprintf(fp, "%f\t%f\n", rsclsgs, rsclvis);
		fclose(fp);

		// ***** modify physical & eddy viscosity accordingly ***** //

		for (int k=1; k<=ms.Nz; k++) {
		for (int i=1; i<=ms.Nx; i++) {
			// boundary strain rate of the coarse velocity field
			const double *sr = vel.ShearStrain(i,j,k);

			double s12 = sr[0], tau12 = tau12sgs(i,j,k) * rsclsgs;
			double s23 = sr[1], tau23 = tau23sgs(i,j,k) * rsclsgs;

			if (k<ms.Nz) nuz(i,j,k) = rsclvis/Re + fmax(.5 * tau12 / s12, 0);
			if (i<ms.Nx) nux(i,j,k) = rsclvis/Re + fmax(.5 * tau23 / s23, 0);
		}}
	}
}


void OffWallSubGridNormal(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double rsclx, double rsclu)
// supply sgs stress 22 on off-wall boundary through viscosity
{
	const Mesh &ms = vis.ms, &ms0 = veldns.ms;

	const Scla &u = vel[1], &u0 = veldns[1];
	const Scla &v = vel[2], &v0 = veldns[2];
	const Scla &w = vel[3], &w0 = veldns[3];

	// viscosity & sgs-stress at cell-center
	Vctr normalsgs(ms);
	const Scla &tau22sgs = normalsgs[2];

	Scla &nuc = vis.GetScl();

	// sgs normal stress filtered from resolved velocity field
	SGS::SubGridNormalStress(normalsgs, veldns, rsclx, rsclu);

	for (int j=1; j<ms.Ny; j+=ms.Ny-2) {

		// the wall normal position on which stresses act
		const double y = ms.yc(j), y0 = WallRscl(y, rsclx);

		// ***** rescale the DNS-filtered sgs stress by Reynolds stress defect ***** //

		double r22dns = 0; // DNS Reynolds stress
		double r22les = 0; // LES resolved Reynolds
		double s22sgs = 0; // mean sgs-stress filtered from DNS

		// calculate reference (DNS) Reynolds shear stress
		int j3 = Interp::BiSearch(y0, ms0.y(), 1, ms0.Ny);
		int j4 = j3 + 1;

		double vm3 = v0.meanxz(j3), y3 = ms0.y(j3);
		double vm4 = v0.meanxz(j4), y4 = ms0.y(j4);

		for (int k=1; k<ms0.Nz; k++) {
		for (int i=1; i<ms0.Nx; i++) {

			r22dns += 1. / (y4-y3) * (
				pow(v0(i,j3,k) - vm3, 2.) * (y4-y0) +
				pow(v0(i,j4,k) - vm4, 2.) * (y0-y3) );
		}}

		r22dns *= rsclu * rsclu
			   * (ms.Nx-1.) / (ms0.Nx-1.) // rescale to match the other two
			   * (ms.Nz-1.) / (ms0.Nz-1.);

		// calculate resolved Reynolds stress and mean sgs-stress
		vm3 = v.meanxz(j3 = j);
		vm4 = v.meanxz(j4 = ms.jpa(j));

		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			s22sgs += tau22sgs(i,j,k);
			r22les += .5 * (
				pow(v(i,j3,k) - vm3, 2.) +
				pow(v(i,j4,k) - vm4, 2.) );
		}}

		const double rsclsgs = fabs((r22dns - r22les) / s22sgs);


		FILE *fp = fopen("WMLOG2.dat", "a");
		if (j==1) fprintf(fp, "%f\t", rsclsgs);
		else      fprintf(fp, "%f\n", rsclsgs);
		fclose(fp);

		// ***** modify physical & eddy viscosity accordingly ***** //

		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {
			// boundary strain rate of the coarse velocity field
			const double *sr = vel.Strainrate(i,j,k);

			double s22 = sr[1], tau22 = tau22sgs(i,j,k) * rsclsgs;

			nuc(i,j,k) = 1./Re + fmax(.5 * tau22 / s22, 0);
		}}
	}
}


void OffWallSubGridDissipation(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double rsclx, double rsclu)
//  modify eddy viscosity near off-wall boundary by sgs dissipation
{
	const Mesh &ms = vis.ms;
	const Mesh &ms0 = veldns.ms;

	Scla &nuc = vis.GetScl();
	Scla &nuz = vis.GetVec(3);

	// DNS strainrate field at cell-centers up to real boundary
	Scla s11(ms0), s22(ms0), s33(ms0);
	Scla s12(ms0), s23(ms0), s13(ms0);

	#pragma omp parallel for
	for (int j=1; j<ms0.Ny; j++) {
	for (int k=1; k<ms0.Nz; k++) {
	for (int i=1; i<ms0.Nx; i++) {

		const double *sr = veldns.Strainrate(i,j,k); // cannot parallelize because sr is static

		s11(i,j,k) = sr[0]; s22(i,j,k) = sr[1]; s33(i,j,k) = sr[2];
		s12(i,j,k) = sr[3]; s23(i,j,k) = sr[4]; s13(i,j,k) = sr[5];
	}}}

	// SGS stress filtered from DNS field
	Vctr shearsgs(ms);
	Vctr normalsgs(ms);

	SGS::SubGridStress(shearsgs, normalsgs, veldns, rsclx, rsclu);

	const Scla &tau12sgs = shearsgs[1], &tau11sgs = normalsgs[1];
	const Scla &tau23sgs = shearsgs[2], &tau22sgs = normalsgs[2];
	const Scla &tau13sgs = shearsgs[3], &tau33sgs = normalsgs[3];

	// re-solve eddy viscosity near off-wall boundary based on sgs dissipation
	for (int j=1; j<ms.Ny; j+=ms.Ny-2) {
	#pragma omp parallel for
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		int id = ms.idx(i,j,k);

		double x = ms.xc(i) * rsclx, dx = ms.dx(i) * rsclx;
		double z = ms.zc(k) * rsclx, dz = ms.dz(k) * rsclx;
		double y = WallRscl(ms.yc(j), rsclx), dy = 0;

		double sgsdsp = rsclu * rsclx * ( // rescale Sij to match inner scale
			tau11sgs[id] * Filter::FilterNodeA(x,y,z,dx,dy,dz,s11) +
			tau22sgs[id] * Filter::FilterNodeA(x,y,z,dx,dy,dz,s22) +
			tau33sgs[id] * Filter::FilterNodeA(x,y,z,dx,dy,dz,s33) + (
			tau12sgs[id] * Filter::FilterNodeA(x,y,z,dx,dy,dz,s12) +
			tau23sgs[id] * Filter::FilterNodeA(x,y,z,dx,dy,dz,s23) +
			tau13sgs[id] * Filter::FilterNodeA(x,y,z,dx,dy,dz,s13) ) * 2 );

		const double *sr = vel.Strainrate(i,j,k); // cannot parallelize because sr is static

		double ss =
			sr[0]*sr[0] + sr[1]*sr[1] + sr[2]*sr[2] +
		  ( sr[3]*sr[3] + sr[4]*sr[4] + sr[5]*sr[5] ) * 2;

		nuc[id] = 1./Re + fmin(fmax(.5*sgsdsp/ss, 0), .1);
	}}}


	// rescale eddy viscosity by Reynolds stress defect
	ProcessBoundaryVis(vis);

	double t12sgs = 0;
	double r12dfc = 0;

	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

		double weight1 = (j>1 ? -1 : 1) * .5 / (ms.Nx-1.) / (ms.Nz-1.);
		double weight2 = (j>1 ? -1 : 1) * .5;

		#pragma omp parallel for reduction(+: t12sgs)
		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			t12sgs += weight1 * 2. * (nuz(i,j,k)-1./Re) * vel.ShearStrain(i,j,k)[0];
		}}

		r12dfc += weight2 * ReynoldsStressDefect(j, vel, veldns, rsclx, rsclu);
	}

	double rscldsp = fmin(fmax(fabs(r12dfc / t12sgs), 1.), 4.);

	FILE *fp = fopen("WMLOG2.dat", "a");
	fprintf(fp, "%f\t%f\n", rscldsp, r12dfc / t12sgs);
	fclose(fp);


	for (int j=1; j<ms.Ny; j+=ms.Ny-2) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		nuc(i,j,k) = 1./Re + rscldsp * (nuc(i,j,k)-1./Re);
	}}}

	ProcessBoundaryVis(vis);
}


void WM::OffWallSGS(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double rsclx, double rsclu)
{
	// OffWallSubGridShear      (vis, vel, veldns, Re, rsclx, rsclu);
	// OffWallSubGridNormal     (vis, vel, veldns, Re, rsclx, rsclu);
	OffWallSubGridDissipation(vis, vel, veldns, Re, rsclx, rsclu);

	// note: when OffWallSubGridDissipation is used in conjuction with OffWallSubGridShear,
	// the boundary processing at the end of OffWallSubGridDissipation should be commented

	// weird bug: OffWallSubGridDissipation must be after OffWallSubGridShear, do not know why...
	// Otherwise NaN may occur, but with a cout of r12dns, r12les, s12sgs, the bug disappears...
}







		// determine rescale factor of sgs dissipation by total turbulent energy discrepency
		// double rscldsp = 2.;

		// double r11dns = 0, r11les = 0;
		// double r22dns = 0, r22les = 0;
		// double r33dns = 0, r33les = 0;

		// const double y0 = WallRscl(ms.yc(j), rsclx);

		// // DNS turbulent energy
		// int j3u = Interp::BiSearch(y0, ms0.yc(), 0, ms0.Ny), j4u = j3u + 1;
		// int j3v = Interp::BiSearch(y0, ms0.y (), 1, ms0.Ny), j4v = j3v + 1;
		// double y3u = ms0.yc(j3u), um3 = u0.meanxz(j3u), wm3 = w0.meanxz(j3u);
		// double y4u = ms0.yc(j4u), um4 = u0.meanxz(j4u), wm4 = w0.meanxz(j4u);
		// double y3v = ms0.y (j3v), vm3 = v0.meanxz(j3v);
		// double y4v = ms0.y (j4v), vm4 = v0.meanxz(j4v);

		// for (int k=1; k<ms0.Nz; k++) {
		// for (int i=1; i<ms0.Nx; i++) {

		// 	r11dns += 1. / (y4u-y3u) * (
		// 		pow(u0(i,j3u,k) - um3, 2.) * (y4u-y0) +
		// 		pow(u0(i,j4u,k) - um4, 2.) * (y0-y3u) );

		// 	r33dns += 1. / (y4u-y3u) * (
		// 		pow(w0(i,j3u,k) - wm3, 2.) * (y4u-y0) +
		// 		pow(w0(i,j4u,k) - wm4, 2.) * (y0-y3u) );

		// 	r22dns += 1. / (y4v-y3v) * (
		// 		pow(v0(i,j3v,k) - vm3, 2.) * (y4v-y0) +
		// 		pow(v0(i,j4v,k) - vm4, 2.) * (y0-y3v) );
		// }}

		// r11dns *= rsclu * rsclu * ((ms.Nx-1.) * (ms.Nz-1.)) / ((ms0.Nx-1.) * (ms0.Nz-1.));
		// r22dns *= rsclu * rsclu * ((ms.Nx-1.) * (ms.Nz-1.)) / ((ms0.Nx-1.) * (ms0.Nz-1.));
		// r33dns *= rsclu * rsclu * ((ms.Nx-1.) * (ms.Nz-1.)) / ((ms0.Nx-1.) * (ms0.Nz-1.));

		// // LES turbulent energy
		// double um = u.meanxz(j);
		// double wm = w.meanxz(j);
		// vm3 = v.meanxz(j3v = j);
		// vm4 = v.meanxz(j4v = ms.jpa(j));

		// for (int k=1; k<ms.Nz; k++) {
		// for (int i=1; i<ms.Nx; i++) {

		// 	r11les += pow(u(i,j,k) - um, 2.);
		// 	r33les += pow(w(i,j,k) - wm, 2.);
		// 	r22les += .5 * (
		// 		pow(v(i,j3v,k) - vm3, 2.) +
		// 		pow(v(i,j4v,k) - vm4, 2.) );
		// }}

		// // rescaling factor
		// rscldsp = (r11dns+r22dns+r33dns) / (r11les+r22les+r33les);

		// FILE *fp = fopen("WMLOG2.dat", "a");
		// if (j==1) fprintf(fp, "%f\t", rscldsp);
		// else      fprintf(fp, "%f\n", rscldsp);
		// fclose(fp);



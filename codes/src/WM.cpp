#include "WM.h"
#include "SGS.h"
#include "Interp.h"

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

void WM::OffWallSGS(
	Flow &vis, const Vctr &vel, const Vctr &veldns,
	double Re, double rsclx, double rsclu)
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
		const double y = ms.y(j);
		const double y0 = y < 1 ? y * rsclx : 2 - (2-y) * rsclx;

		// ***** rescale the DNS-filtered sgs stress by Reynolds stress defect ***** //

		double r12dns = 0; // DNS Reynolds stress
		double r12les = 0; // LES resolved Reynolds
		double s12sgs = 0; // mean sgs-stress filtered from DNS

		// calculate reference (DNS) Reynolds shear stress
		int j3u = Interp::BiSearch(y0, ms0.yc(), 0, ms0.Ny), j4u = j3u+1;
		int j3v = Interp::BiSearch(y0, ms0.y(),  1, ms0.Ny), j4v = j3v+1;

		double um3 = u0.meanxz(j3u), y3u = ms0.yc(j3u);
		double um4 = u0.meanxz(j4u), y4u = ms0.yc(j4u);
		double vm3 = v0.meanxz(j3v), y3v = ms0.y(j3v);
		double vm4 = v0.meanxz(j4v), y4v = ms0.y(j4v);

		double um = (um3 * (y4u-y0) + um4 * (y0-y3u)) / (y4u-y3u);
		double vm = (vm3 * (y4v-y0) + vm4 * (y0-y3v)) / (y4v-y3v);

		#pragma omp parallel for reduction(+: r12dns)
		for (int k=1; k<ms0.Nz; k++) { double z = ms0.zc(k);
		for (int i=1; i<ms0.Nx; i++) { double x = ms0.x(i);
			// calculate uv on z-edges and interpolate to the desired y position
			r12dns += (Interp::InterpNodeU(x,y0,z,u0) - um)
			        * (Interp::InterpNodeV(x,y0,z,v0) - vm);
		}}

		r12dns *= rsclu * rsclu
			   * (ms.Nx-1.) / (ms0.Nx-1.) // rescale to match the other two
			   * (ms.Nz-1.) / (ms0.Nz-1.);

		// calculate resolved Reynolds stress and mean sgs-stress
		um3 = u.meanxz(j3u = ms.jma(j)); y3u = ms.yc(j3u);
		um4 = u.meanxz(j4u = j);         y4u = ms.yc(j4u);
		vm  = v.meanxz(j);
		um  = (um3 * (y4u-y) + um4 * (y-y3u)) / (y4u-y3u);

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

		const double rsclsgs = fabs((r12dns - r12les) / s12sgs);
		const double rsclvis = fabs((y4u-y3u) / (1-fabs(y-1)) / log((1-fabs(y4u-1))/(1-fabs(y3u-1))));


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

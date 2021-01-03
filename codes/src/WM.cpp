#include "WM.h"
#include "SGS.h"
#include "Interp.h"
#include "Filter.h"
#include "Bcond.h"
#include "PIO.h"
#include "Statis.h"

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
		nuz(i,j,k) = fmax((j==1 ? 1 : -1) * .5*tau12/s12, 0);
	}}}
}


void WM::LogLawWallShear(Flow &vis, const Vctr &vel, double Ret)
// determine the local wall shear by log-law, as proposed by Piomelli et al. (1989)
{
	const Mesh &ms = vis.ms;
	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	Scla &nuz = vis.GetVec(3);
	Scla &nux = vis.GetVec(1);

	const double kappa = .41;
	const double B     = 5.3;
	const double theta = 13. * PI/180;

	double y = ms.yc(2);
	double law = 1. / (log(y * Ret) / kappa + B);
	double dlt = y / tan(theta);

	#pragma omp parallel for collapse(3)
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double yyy = j==1 ? y : 2-y;
		double sgn = j==1 ? 1 : -1;

		double tau12 = sgn * law * Interp::InterpNodeU(dlt+ms.x(i),yyy,ms.zc(k),u);
		double tau23 = sgn * law * Interp::InterpNodeW(dlt+ms.xc(i),yyy,ms.z(k),w);

		const double *sr = vel.ShearStrain(i,j,k);

		nuz(i,j,k) = fmin(fmax(.5 * tau12 / sr[0], -.1), .1);
		nux(i,j,k) = fmin(fmax(.5 * tau23 / sr[1], -.1), .1);
	}}}
}


void WM::OffWallSubGridUniform(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu)
{
	const Mesh &ms = vis.ms;

	Scla &nux = vis.GetVec(1);
	Scla &nuz = vis.GetVec(3);

	// Reynolds stress defect
	double r12dfc = veldns.ms.Ly < 2 ? \
		// HALFMFU
		ReynoldsStressDefect(1, vel, veldns, Re, Ret, rsclx, rsclu) \
		// full channel
		: .5 * (
		ReynoldsStressDefect(1,     vel, veldns, Re, Ret, rsclx, rsclu) -
		ReynoldsStressDefect(ms.Ny, vel, veldns, Re, Ret, rsclx, rsclu) );

	// r12dfc = fmin(r12dfc, 0); // positive mean SGS stress is forbiden

	FILE* fp = fopen("r12dfc.dat", "a");
	fprintf(fp, "%.6e\n", r12dfc);
	fclose(fp);

	// rescale kinematic viscosity
	const double y = ms.y(1), y3 = ms.yc(0), y4 = ms.yc(1);
	const double rsclvis = 1.; //fabs((y4-y3) / (1-fabs(y-1)) / log((1-fabs(y4-1))/(1-fabs(y3-1))));

	// ***** modify kinematic & eddy viscosity accordingly ***** //
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

		#pragma omp parallel for // ShearStrain() returns threadprivate pointer, so it is parallel-safe
		for (int k=1; k<=ms.Nz; k++) {
		for (int i=1; i<=ms.Nx; i++) {
			// boundary strain rate of the coarse velocity field
			const double *sr = vel.ShearStrain(i,j,k);

			double s12 = sr[0], tau12 = r12dfc * (j<=1 ? -1 : 1);
			double s23 = sr[1], tau23 = 0;

			// it is necessary to allow negative eddy viscosity on boundary for correct tau
			if (k<ms.Nz) nuz(i,j,k) = rsclvis/Re + fmin(fmax(.5 * tau12 / s12, -.1), .1);
			if (i<ms.Nx) nux(i,j,k) = rsclvis/Re + fmin(fmax(.5 * tau23 / s23, -.1), .1);
		}}
	}
}


void WM::OffWallSubGridShear(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu)
// supply sgs stress 12 & 23 on off-wall boundary through viscosity
{
	const Mesh &ms = vis.ms;

	// viscosity & sgs-stress on edges aimed for
	Vctr shearsgs(ms);
	const Scla &tau23sgs = shearsgs[2];
	const Scla &tau12sgs = shearsgs[1];

	Scla &nux = vis.GetVec(1);
	Scla &nuz = vis.GetVec(3);

	// sgs shear stress filtered from resolved velocity field
	SGS::SubGridShearStress(shearsgs, veldns, rsclx, rsclu);

	PIO::PredictBoundarySGS(shearsgs, vel, Ret);

	// ***** rescale the DNS-filtered sgs stress by Reynolds stress defect ***** //
	double r12dfc = 0; // Reynolds stress defect
	double t12sgs = 0; // mean SGS stress (component 12)

	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
		// average for top & bottom real boundaries
		double weight0 = (j>1 ? -1 : 1) * .5 + .5 * (veldns.ms.Ly < 2); // handle HALFMFU
		double weight1 = (j>1 ? -1 : 1) * .5 / (ms.Nx-1.) / (ms.Nz-1.);

		r12dfc += weight0 * ReynoldsStressDefect(j, vel, veldns, Re, Ret, rsclx, rsclu);

		for (int k=1; k<ms.Nz; k++)
		for (int i=1; i<ms.Nx; i++)
			t12sgs += weight1 * tau12sgs(i,j,k);
	}

	const double rsclsgs = fmin(fabs(r12dfc / t12sgs), 2.);

	// ***** rescale kinematic viscosity to account for low-order differencing error ***** //
	const double y = ms.y(1), y3 = ms.yc(0), y4 = ms.yc(1);
	const double rsclvis = 1.; //fabs((y4-y3) / (1-fabs(y-1)) / log((1-fabs(y4-1))/(1-fabs(y3-1))));

	FILE *fp = fopen("WMLOG.dat", "a");
	fprintf(fp, "%f\t%f\n", rsclsgs, rsclvis);
	fclose(fp);

	// ***** modify kinematic & eddy viscosity accordingly ***** //
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

		#pragma omp parallel for // ShearStrain() returns threadprivate pointer, so it is parallel-safe
		for (int k=1; k<=ms.Nz; k++) {
		for (int i=1; i<=ms.Nx; i++) {
			// boundary strain rate of the coarse velocity field
			const double *sr = vel.ShearStrain(i,j,k);

			// double s12 = sr[0], tau12 = tau12sgs(i,j,k) * rsclsgs;
			// double s23 = sr[1], tau23 = tau23sgs(i,j,k) * rsclsgs;

			double s12 = sr[0], tau12 = tau12sgs(i,j,k) + t12sgs * (rsclsgs - 1) * (j==1 ? 1 : -1);
			double s23 = sr[1], tau23 = tau23sgs(i,j,k);

			// it is necessary to allow negative eddy viscosity on boundary for correct tau
			if (k<ms.Nz) nuz(i,j,k) = rsclvis/Re + fmin(fmax(.5 * tau12 / s12, -.1), .1);
			if (i<ms.Nx) nux(i,j,k) = rsclvis/Re + fmin(fmax(.5 * tau23 / s23, -.1), .1);
		}}
	}
}


void WM::OffWallSubGridDissipation(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu)
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
		double y = Filter::WallRscl(ms.yc(j), rsclx), dy = 0;

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
	{
		Scla &nuc = vis.GetScl();
		Bcond::SetBoundaryY(nuc, 1); // homogeneous Neumann
		Bcond::SetBoundaryX(nuc, 3); // periodic
		Bcond::SetBoundaryZ(nuc, 3); // periodic
		vis.CellCenter2Edge();
	}

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

		r12dfc += weight2 * ReynoldsStressDefect(j, vel, veldns, Re, Ret, rsclx, rsclu);
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

	{
		Scla &nuc = vis.GetScl();
		Bcond::SetBoundaryY(nuc, 1); // homogeneous Neumann
		Bcond::SetBoundaryX(nuc, 3); // periodic
		Bcond::SetBoundaryZ(nuc, 3); // periodic
		vis.CellCenter2Edge();
	}

	// note: when OffWallSubGridDissipation is used in conjuction with OffWallSubGridShear,
	// the boundary processing at the end of OffWallSubGridDissipation should be commented

	// weird bug: OffWallSubGridDissipation must be after OffWallSubGridShear, do not know why...
	// Otherwise NaN may occur, but with a cout of r12dns, r12les, s12sgs, the bug disappears...
}



void WM::debug_ShowBoundaryShear(const Vctr &vel, const Flow &vis)
{
	const Mesh &ms = vel.ms;

	double shearrey = 0;
	double shearvis = 0;

	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
		shearrey += (j==1 ? -1 : 1) * .5 * Statis::ReynoldsStress(j, vel);
		shearvis += (j==1 ? 1 : -1) * .5 * Statis::MeanVisShearStresses (j, vel, vis)[0];
	}

	cout << shearvis + shearrey << '\t' << shearrey << '\t' << shearvis << endl << endl;
}



double WM::ReynoldsStressDefect(int j, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu)
// calculate Reynolds stress R12 - R12^r at Z-edge of layer j
{
	const Mesh &ms = vel.ms;
	const Mesh &ms0= veldns.ms;

	// the wall normal position on which stresses act
	const double y = ms.y(j);
	const double y0 = Filter::WallRscl(y, rsclx); // position in DNS field with y^+ matched

	// calculate reference (DNS) Reynolds shear stress
	int j3 = Interp::BiSearch(y0, ms0.y(), 1, ms0.Ny);
	int j4 = j3 + 1;

	double y3 = ms0.y(j3);
	double y4 = ms0.y(j4);

	double r12dns = (
		Statis::ReynoldsStress(j3, veldns) * (y4-y0)
	  + Statis::ReynoldsStress(j4, veldns) * (y0-y3)) / (y4-y3);
	
	r12dns *= rsclu * rsclu;
	r12dns += pow(Ret/Re, 2.) * (1 - fabs(y-1.)) * (1 - rsclx);

	// LES resolved Reynolds
	double r12les = Statis::ReynoldsStress(j, vel);

	return r12dns - r12les;
}









// a test implementation that has been proved of no use
void OffWallSubGridShear_from_inner(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu)
// rescale SGS stress 12 & 23 on off-wall boundary through viscosity
{
	const Mesh &ms = vis.ms;
	Scla &nux = vis.GetVec(1);
	Scla &nuz = vis.GetVec(3);

	// ***** rescale boundary SGS stress by Reynolds stress defect ***** //
	double r12dfc = 0;
	double t12sgs = 0;

	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

		double weight0 = (j>1 ? -1 : 1) * .5;
		double weight1 = (j>1 ? -1 : 1) * .5 / (ms.Nx-1.) / (ms.Nz-1.);

		r12dfc += weight0 * WM::ReynoldsStressDefect(j, vel, veldns, Re, Ret, rsclx, rsclu);

		#pragma omp parallel for reduction(+: t12sgs)
		for (int k=1; k<ms.Nz; k++)
		for (int i=1; i<ms.Nx; i++)
			t12sgs += weight1 * 2. * (nuz(i,j,k)-1./Re) * vel.ShearStrain(i,j,k)[0];
	}

	const double rsclsgs = fmin(fabs(r12dfc / t12sgs), 4.);

	// ***** rescale kinematic viscosity to account for low-order differencing error ***** //
	const double y = ms.y(1), y3 = ms.yc(0), y4 = ms.yc(1);
	const double rsclvis = fabs((y4-y3) / (1-fabs(y-1)) / log((1-fabs(y4-1))/(1-fabs(y3-1))));

	FILE *fp = fopen("WMLOG.dat", "a");
	fprintf(fp, "%f\t%f\n", rsclsgs, rsclvis);
	fclose(fp);

	// ***** modify kinematic & eddy viscosity accordingly ***** //
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

		#pragma omp parallel for
		for (int k=1; k<=ms.Nz; k++) {
		for (int i=1; i<=ms.Nx; i++) {

			nuz(i,j,k) = rsclvis/Re + rsclsgs * (nuz(i,j,k)-1./Re);
			nux(i,j,k) = rsclvis/Re + rsclsgs * (nux(i,j,k)-1./Re);
		}}
	}
}



		// // debug: record s12 & tau12 for JPDF
		// if (j == 1) {
		// 	static int cnt = 0;
		// 	if (cnt % 100 == 0) {
		// 		char str[128];
		// 		sprintf(str, "records/tau_record%d.txt", cnt);
		// 		FILE *fp = fopen(str, "w");
		// 		for (int k=1; k<=ms.Nz; k++) {
		// 		for (int i=1; i<=ms.Nx; i++) {
		// 			fprintf(fp, "%.18e\t%.18e\n", tau12sgs(i,j,k), vel.ShearStrain(i,j,k)[0]);
		// 		}}
		// 		fclose(fp);
		// 	}
		// 	cnt ++;
		// }





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








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
		nuz(i,j,k) = fmax((j==1 ? 1 : -1) * .5*tau12/s12, 0);
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

	double um3 = 0, vm3 = 0, y3 = ms0.yc(j3);
	double um4 = 0, vm4 = 0, y4 = ms0.yc(j4);

	#pragma omp parallel for reduction(+: r12dns,um3,vm3,um4,vm4)
	for (int k=1; k<ms0.Nz; k++) {
	for (int i=1; i<ms0.Nx; i++) {

		double ua3 = .5 * (u0(i,j3,k) + u0(ms0.ipa(i),j3,k));
		double va3 = .5 * (v0(i,j3,k) + v0(i,ms0.jpa(j3),k));
		double ua4 = .5 * (u0(i,j4,k) + u0(ms0.ipa(i),j4,k));
		double va4 = .5 * (v0(i,j4,k) + v0(i,ms0.jpa(j4),k));

		um3 += ua3; vm3 += va3;
		um4 += ua4; vm4 += va4;
		r12dns += (ua3*va3 * (y4-y0) + ua4*va4 * (y0-y3)) / (y4-y3);
	}}

	double weight = 1. / (ms0.Nx-1.) / (ms0.Nz-1.);

	r12dns -= weight * (um3*vm3 * (y4-y0) + um4*vm4 * (y0-y3)) / (y4-y3);
	r12dns *= weight * rsclu * rsclu;

	// calculate resolved Reynolds stress
	double um = 0;
	double vm = 0;

	#pragma omp parallel for reduction(+: r12les,um,vm)
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double ua = .5/ms.hy(j) * (u(i,j,k) * ms.dy(j-1) + u(i,ms.jma(j),k) * ms.dy(j));
		double va = .5/ms.hx(i) * (v(i,j,k) * ms.dx(i-1) + v(ms.ima(i),j,k) * ms.dx(i));

		um += ua;
		vm += va;
		r12les += ua * va;
	}}

	weight = 1. / (ms.Nx-1.) / (ms.Nz-1.);

	r12les -= weight * um * vm;
	r12les *= weight;

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
	const Mesh &ms = vis.ms;

	// viscosity & sgs-stress on edges aimed for
	Vctr shearsgs(ms);
	const Scla &tau23sgs = shearsgs[2];
	const Scla &tau12sgs = shearsgs[1];

	Scla &nux = vis.GetVec(1);
	Scla &nuz = vis.GetVec(3);

	// sgs shear stress filtered from resolved velocity field
	SGS::SubGridShearStress(shearsgs, veldns, rsclx, rsclu);

	// ***** rescale the DNS-filtered sgs stress by Reynolds stress defect ***** //
	double r12dfc = 0; // Reynolds stress defect
	double t12sgs = 0; // mean SGS stress (component 12)

	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
		// average for top & bottom real boundaries
		double weight0 = (j>1 ? -1 : 1) * .5;
		double weight1 = (j>1 ? -1 : 1) * .5 / (ms.Nx-1.) / (ms.Nz-1.);

		r12dfc += weight0 * ReynoldsStressDefect(j, vel, veldns, rsclx, rsclu);

		for (int k=1; k<ms.Nz; k++)
		for (int i=1; i<ms.Nx; i++)
			t12sgs += weight1 * tau12sgs(i,j,k);
	}

	const double rsclsgs = fmin(fmax( fabs(r12dfc / t12sgs), 1.), 2.);

	// ***** rescale kinematic viscosity to account for low-order differencing error ***** //
	const double y = ms.y(1), y3 = ms.yc(0), y4 = ms.yc(1);
	const double rsclvis = fabs((y4-y3) / (1-fabs(y-1)) / log((1-fabs(y4-1))/(1-fabs(y3-1))));

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

			double s12 = sr[0], tau12 = tau12sgs(i,j,k) * rsclsgs;
			double s23 = sr[1], tau23 = tau23sgs(i,j,k) * rsclsgs;

			// it is necessary to allow negative eddy viscosity on boundary for correct tau
			if (k<ms.Nz) nuz(i,j,k) = rsclvis/Re + fmin(fmax(.5 * tau12 / s12, -.1), .1);
			if (i<ms.Nx) nux(i,j,k) = rsclvis/Re + fmin(fmax(.5 * tau23 / s23, -.1), .1);
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
	OffWallSubGridShear      (vis, vel, veldns, Re, rsclx, rsclu);
	// OffWallSubGridNormal     (vis, vel, veldns, Re, rsclx, rsclu);
	// OffWallSubGridDissipation(vis, vel, veldns, Re, rsclx, rsclu);

	// note: when OffWallSubGridDissipation is used in conjuction with OffWallSubGridShear,
	// the boundary processing at the end of OffWallSubGridDissipation should be commented

	// weird bug: OffWallSubGridDissipation must be after OffWallSubGridShear, do not know why...
	// Otherwise NaN may occur, but with a cout of r12dns, r12les, s12sgs, the bug disappears...
}





void WM::debug_ShowBoundaryShear(const Vctr &vel, const Flow &vis)
{
	const Mesh &ms = vel.ms;
	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &nuz = vis.SeeVec(3);

	double shearrey = 0;
	double shearvis = 0;

	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

		double um = 0;
		double vm = 0;
		double r12 = 0;
		double t12 = 0;

		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			double ua = .5/ms.hy(j) * (u(i,j,k) * ms.dy(j-1) + u(i,ms.jma(j),k) * ms.dy(j));
			double va = .5/ms.hx(i) * (v(i,j,k) * ms.dx(i-1) + v(ms.ima(i),j,k) * ms.dx(i));

			um += ua;
			vm += va;
			r12 += ua * va;
			t12 += 2 * nuz(i,j,k) * vel.ShearStrain(i,j,k)[0];
		}}

		double weight = 1. / (ms.Nx-1.) / (ms.Nz-1.);

		r12 -= weight * um * vm;
		r12 *= weight;
		t12 *= weight;

		shearrey += (j==1 ? 1 : -1) * .5 * ( - r12);
		shearvis += (j==1 ? 1 : -1) * .5 * t12;
	}

	cout << shearvis + shearrey << '\t' << shearrey << '\t' << shearvis << endl << endl;
}










// a test implementation that has been proved of no use
void OffWallSubGridShear_from_inner(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double rsclx, double rsclu)
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

		r12dfc += weight0 * ReynoldsStressDefect(j, vel, veldns, rsclx, rsclu);

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








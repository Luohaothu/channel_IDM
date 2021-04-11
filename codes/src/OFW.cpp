#include "OFW.h"
#include "Interp.h"
#include "Filter.h"
#include "Statis.h"
#include "PIO.h"
#include "SGS.h"
#include "WM.h"

using namespace std;


void OFW::OffWallSubGridUniform(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu)
{
	const Mesh &ms = vis.ms;
	const Mesh &ms0 = veldns.ms;

	double r23ref = 0;
	double r12ref = Statis::ReynoldsStress(
		Interp::BiSearch(ms.y(1)*rsclx, ms0.y(), 1, ms0.Ny), veldns
		) * rsclu * rsclu + pow(Ret/Re, 2.) * ms.y(1) * (1 - rsclx); // in case MFU is at a different Re than LES

	WM::UniformReyShear(vis, vel, r12ref, r23ref, Re);
}

void OFW::OffWallSubGridShear(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu)
// supply sgs stress 12 & 23 on off-wall boundary through viscosity
{
	const Mesh &ms = vis.ms;
	const Mesh &ms0 = veldns.ms;

	Vctr shearsgs(ms); // sgs-stress on edges
	Scla &tau12 = shearsgs[1];
	Scla &tau23 = shearsgs[2];

	// predicted sgs shear stress (filter + modulation)
	SGS::SubGridShearStress(shearsgs, veldns, rsclx, rsclu);
	PIO::PredictBoundarySGS(shearsgs, vel, Ret);

	double r23ref = 0;
	double r12ref = Statis::ReynoldsStress(
		Interp::BiSearch(ms.y(1)*rsclx, ms0.y(), 1, ms0.Ny), veldns
		) * rsclu * rsclu + pow(Ret/Re, 2.) * ms.y(1) * (1 - rsclx); // in case MFU is at a different Re than LES

	const vector<double> rs3 = Statis::ReynoldsShearStresses(1,     vel);
	const vector<double> rs4 = Statis::ReynoldsShearStresses(ms.Ny, vel);

	double tau12m = .5 * (rs3[0]-rs4[0]);
	double tau23m = .5 * (rs3[1]-rs4[1]);
	double weight = .5 / (ms.Nx-1.) / (ms.Nz-1.);

	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		tau12m -= weight * (tau12(i,1,k) - tau12(i,ms.Ny,k));
		tau23m -= weight * (tau23(i,1,k) - tau23(i,ms.Ny,k));
	}}

	tau12.AddLyr(tau12m - r12ref, 1).AddLyr(r12ref - tau12m, ms.Ny);
	tau23.AddLyr(tau23m - r23ref, 1).AddLyr(r23ref - tau23m, ms.Ny);

	WM::SubGridShear(vis, vel, shearsgs, Re);
}

void OFW::OffWallSubGridDissipation(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu)
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
		#pragma omp parallel
		vis.CellCenter2Edges();
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
		#pragma omp parallel
		vis.CellCenter2Edges();
	}

	// note: when OffWallSubGridDissipation is used in conjuction with OffWallSubGridShear,
	// the boundary processing at the end of OffWallSubGridDissipation should be commented

	// weird bug: OffWallSubGridDissipation must be after OffWallSubGridShear, do not know why...
	// Otherwise NaN may occur, but with a cout of r12dns, r12les, s12sgs, the bug disappears...
}


void OFW::OffWallVelo(Boundaries &bc, Boundaries &sbc, Vctr &veltmp,
	const Vctr &vel, const Vctr &velmfu, double Ret, double rsclx, double rsclu)
{
	PIO::PredictBoundary(veltmp, vel, velmfu, Ret, rsclx, rsclu);
	Bcond::ChannelDirichlet(bc, sbc, vel.ms, veltmp);
}


double OFW::ReynoldsStressDefect(int j, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu)
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



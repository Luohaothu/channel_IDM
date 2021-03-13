#include "WM.h"
#include "SGS.h"
#include "Interp.h"
#include "Filter.h"
#include "Bcond.h"
#include "Statis.h"

using namespace std;


void WM::UniformVisShear(Flow &vis, const Vctr &vel, double tau12, double tau23)
// modify viscosity on x&z-edges to ensure uniform (kinematic + eddy) viscous shear shress on every boundary point
{
	const Mesh &ms = vis.ms;

	Scla &nux = vis.GetVec(1);
	Scla &nuz = vis.GetVec(3);

	#pragma omp parallel for collapse(3)
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
	for (int k=1; k<=ms.Nz; k++) {
	for (int i=1; i<=ms.Nx; i++) {

		double sgn = j==1 ? 1 : -1;
		const double *sr = vel.ShearStrain(i,j,k);

		if (k<ms.Nz) nuz(i,j,k) = fmin(fmax(sgn * .5 * tau12 / sr[0], -.1), .1);
		if (i<ms.Nx) nux(i,j,k) = fmin(fmax(sgn * .5 * tau23 / sr[1], -.1), .1);
	}}}
}

void WM::UniformReyShear(Flow &vis, const Vctr &vel, double r12ref, double r23ref, double Re)
// modify viscosity on x&z-edges to ensure uniform (SGS + resolvable) Reynolds shear stress on every boundary point
{
	const Mesh &ms = vis.ms;

	Scla &nux = vis.GetVec(1);
	Scla &nuz = vis.GetVec(3);

	const vector<double> rs3 = Statis::ReynoldsShearStresses(1,     vel);
	const vector<double> rs4 = Statis::ReynoldsShearStresses(ms.Ny, vel);
	
	double tau12 = .5 * (rs3[0]-rs4[0]) - r12ref; // defect in Reynolds stress
	double tau23 = .5 * (rs3[1]-rs4[1]) - r23ref;

	double rsclvis = 1.; // support some adjustment (such as log) to mean velocity

	#pragma omp parallel for collapse(3)
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
	for (int k=1; k<=ms.Nz; k++) {
	for (int i=1; i<=ms.Nx; i++) {

		double sgn = j==1 ? 1 : -1;
		const double *sr = vel.ShearStrain(i,j,k);

		if (k<ms.Nz) nuz(i,j,k) = rsclvis/Re + fmin(fmax(sgn * .5 * tau12 / sr[0], -.1), .1);
		if (i<ms.Nx) nux(i,j,k) = rsclvis/Re + fmin(fmax(sgn * .5 * tau23 / sr[1], -.1), .1);
	}}}

	FILE* fp = fopen("r12dfc.dat", "a");
	fprintf(fp, "%.6e\n", tau12);
	fclose(fp);
}

void WM::LogLawWallShear(Flow &vis, const Vctr &vel, double Ret)
// determine the local wall shear by log-law, as proposed by Piomelli et al. (1989)
{
	const Mesh &ms = vis.ms;
	const Scla &u = vel[1];
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
	for (int k=1; k<=ms.Nz; k++) {
	for (int i=1; i<=ms.Nx; i++) {

		double yyy = j==1 ? y : 2-y;
		double sgn = j==1 ? 1 : -1;

		double tau12 = sgn * law * Interp::InterpNodeU(dlt+ms.x(i),yyy,ms.zc(k),u);
		double tau23 = sgn * law * Interp::InterpNodeW(dlt+ms.xc(i),yyy,ms.z(k),w);

		const double *sr = vel.ShearStrain(i,j,k);

		if (k<ms.Nz) nuz(i,j,k) = fmin(fmax(.5 * tau12 / sr[0], -.1), .1);
		if (i<ms.Nx) nux(i,j,k) = fmin(fmax(.5 * tau23 / sr[1], -.1), .1);
	}}}
}

void WM::SubGridShear(Flow &vis, const Vctr &vel, const Vctr &shearsgs, double Re)
// supply sgs stress 12 & 23 on off-wall boundary through viscosity
{
	const Mesh &ms = vis.ms;
	const Scla &tau12 = shearsgs[1]; // sgs-stress on edges aimed for
	const Scla &tau23 = shearsgs[2];

	Scla &nux = vis.GetVec(1);
	Scla &nuz = vis.GetVec(3);

	double rsclvis = 1.; // support some adjustment (such as log) to mean velocity

	#pragma omp parallel for collapse(3)
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
	for (int k=1; k<=ms.Nz; k++) {
	for (int i=1; i<=ms.Nx; i++) {

		const double *sr = vel.ShearStrain(i,j,k);

		// it is necessary to allow negative eddy viscosity on boundary for correct tau
		if (k<ms.Nz) nuz(i,j,k) = rsclvis/Re + fmin(fmax(.5 * tau12(i,j,k) / sr[0], -.1), .1);
		if (i<ms.Nx) nux(i,j,k) = rsclvis/Re + fmin(fmax(.5 * tau23(i,j,k) / sr[1], -.1), .1);
	}}}
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









// // a test implementation that has been proved of no use
// void OffWallSubGridShear_from_inner(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu)
// // rescale SGS stress 12 & 23 on off-wall boundary through viscosity
// {
// 	const Mesh &ms = vis.ms;
// 	Scla &nux = vis.GetVec(1);
// 	Scla &nuz = vis.GetVec(3);

// 	// ***** rescale boundary SGS stress by Reynolds stress defect ***** //
// 	double r12dfc = 0;
// 	double t12sgs = 0;

// 	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

// 		double weight0 = (j>1 ? -1 : 1) * .5;
// 		double weight1 = (j>1 ? -1 : 1) * .5 / (ms.Nx-1.) / (ms.Nz-1.);

// 		r12dfc += weight0 * WM::ReynoldsStressDefect(j, vel, veldns, Re, Ret, rsclx, rsclu);

// 		#pragma omp parallel for reduction(+: t12sgs)
// 		for (int k=1; k<ms.Nz; k++)
// 		for (int i=1; i<ms.Nx; i++)
// 			t12sgs += weight1 * 2. * (nuz(i,j,k)-1./Re) * vel.ShearStrain(i,j,k)[0];
// 	}

// 	const double rsclsgs = fmin(fabs(r12dfc / t12sgs), 4.);

// 	// ***** rescale kinematic viscosity to account for low-order differencing error ***** //
// 	const double y = ms.y(1), y3 = ms.yc(0), y4 = ms.yc(1);
// 	const double rsclvis = fabs((y4-y3) / (1-fabs(y-1)) / log((1-fabs(y4-1))/(1-fabs(y3-1))));

// 	FILE *fp = fopen("WMLOG.dat", "a");
// 	fprintf(fp, "%f\t%f\n", rsclsgs, rsclvis);
// 	fclose(fp);

// 	// ***** modify kinematic & eddy viscosity accordingly ***** //
// 	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

// 		#pragma omp parallel for
// 		for (int k=1; k<=ms.Nz; k++) {
// 		for (int i=1; i<=ms.Nx; i++) {

// 			nuz(i,j,k) = rsclvis/Re + rsclsgs * (nuz(i,j,k)-1./Re);
// 			nux(i,j,k) = rsclvis/Re + rsclsgs * (nux(i,j,k)-1./Re);
// 		}}
// 	}
// }



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








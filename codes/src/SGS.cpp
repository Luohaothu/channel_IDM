#include "SGS.h"
#include "Bcond.h"
#include "Filter.h"

using namespace std;


// SGS::SGS(const Mesh &ms):
// s11(ms), s22(ms), s33(ms), s12(ms), s23(ms), s13(ms),
// m11(ms), m22(ms), m33(ms), m12(ms), m23(ms), m13(ms),
// l11(ms), l22(ms), l33(ms), l12(ms), l23(ms), l13(ms)
// {}

// note: MUST NOT filter and assign values to an array in the same transverse


void SGS::Smargorinsky(Scla &nut, const Vctr &vel, double Re, double Cs)
{
	const Mesh &ms = nut.ms;
	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	// calculate delta_nu. NOTE: this only applys when boundary is on the wall
	double tauw = 0;

	for (int j=1; j<=ms.Ny; j+= ms.Ny-1) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double s12 = .5 * (
			1./ms.hx(i) * (v(i,j,k) - v(ms.ima(i),j,k)) +
			1./ms.hy(j) * (u(i,j,k) - u(i,ms.jma(j),k)) );

		double weight = (j==1 ? .5 : -.5) / (ms.Nx-1) / (ms.Nz-1);

		tauw += 2./Re * s12 * weight;
	}}}

	double Ret = Re * sqrt(tauw); // Ly is taken for 2.0 as default

	// calculate eddy viscosity at cell centers
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double dlt = pow(ms.dx(i)*ms.dy(j)*ms.dz(k), 1./3); // filter size
		double dmp = 1 - exp( (fabs(1-ms.yc(j)) - 1) * Ret / 25. ); // Van Driest damping: 1-exp(-y^+/A^+), with A^+ = 25

		const double *sr = vel.Strainrate(i,j,k);

		double sra = sqrt(
			2. * ( sr[0]*sr[0] + sr[1]*sr[1] + sr[2]*sr[2] +
			2. * ( sr[3]*sr[3] + sr[4]*sr[4] + sr[5]*sr[5] )));

		nut(i,j,k) = pow(Cs * dlt * dmp, 2.) * sra;
	}}}

	// set boundary eddy viscosity
	Bcond::SetBoundaryY(nut, 1); // homogeneous Neumann
	Bcond::SetBoundaryX(nut, 3); // periodic
	Bcond::SetBoundaryZ(nut, 3); // periodic
}

void SGS::DynamicSmarg(Scla &nut, const Vctr &vel)
{
	const Mesh &ms = nut.ms;
	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	Scla s11(ms), s22(ms), s33(ms), s12(ms), s23(ms), s13(ms);
	Scla m11(ms), m22(ms), m33(ms), m12(ms), m23(ms), m13(ms);
	Scla l11(ms), l22(ms), l33(ms), l12(ms), l23(ms), l13(ms);
	Scla &uf=l11, &vf=l22, &wf=l33, &uc=l12, &vc=l23, &wc=l13;

	/***** solve Mij *****/

#pragma omp parallel
{
	// Sij, |Sij| and |Sij|*Sij stored in Lij
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		const double *sr = vel.Strainrate(i,j,k); // returned value in Strainrate is static and may not be parallised

		double sra = sqrt(
			2. * ( sr[0]*sr[0] + sr[1]*sr[1] + sr[2]*sr[2] +
			2. * ( sr[3]*sr[3] + sr[4]*sr[4] + sr[5]*sr[5] )));

		// Sij
		s11(i,j,k) = sr[0]; s22(i,j,k) = sr[1]; s33(i,j,k) = sr[2];
		s12(i,j,k) = sr[3]; s23(i,j,k) = sr[4]; s13(i,j,k) = sr[5];

		// |Sij|*Sij, stored in Lij
		l11(i,j,k) = sra * sr[0]; l12(i,j,k) = sra * sr[3];
		l22(i,j,k) = sra * sr[1]; l23(i,j,k) = sra * sr[4];
		l33(i,j,k) = sra * sr[2]; l13(i,j,k) = sra * sr[5];

		// |Sij|
		nut(i,j,k) = sra;
	}}}

	// Mij
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		// filter & test filter sizes
		double dlt1_sq = pow(ms.dx(i)*ms.dy(j)*ms.dz(k), 2./3);
		double dlt2_sq = pow(2., 4./3) * dlt1_sq;

		// F(Sij)
		double fsr[6];

		fsr[0] = Filter::TestFilter(i,j,k,s11);
		fsr[1] = Filter::TestFilter(i,j,k,s22);
		fsr[2] = Filter::TestFilter(i,j,k,s33);
		fsr[3] = Filter::TestFilter(i,j,k,s12);
		fsr[4] = Filter::TestFilter(i,j,k,s23);
		fsr[5] = Filter::TestFilter(i,j,k,s13);

		// |F(Sij)|
		double fsra = sqrt(
			2. * ( fsr[0]*fsr[0] + fsr[1]*fsr[1] + fsr[2]*fsr[2] +
			2. * ( fsr[3]*fsr[3] + fsr[4]*fsr[4] + fsr[5]*fsr[5] )));

		// Mij = dlt_2^2 * |F(Sij)| * F(Sij) - dlt_1^2 * F(|Sij|*Sij)
		m11(i,j,k) = dlt2_sq * fsra * fsr[0] - dlt1_sq * Filter::TestFilter(i,j,k,l11);
		m22(i,j,k) = dlt2_sq * fsra * fsr[1] - dlt1_sq * Filter::TestFilter(i,j,k,l22);
		m33(i,j,k) = dlt2_sq * fsra * fsr[2] - dlt1_sq * Filter::TestFilter(i,j,k,l33);
		m12(i,j,k) = dlt2_sq * fsra * fsr[3] - dlt1_sq * Filter::TestFilter(i,j,k,l12);
		m23(i,j,k) = dlt2_sq * fsra * fsr[4] - dlt1_sq * Filter::TestFilter(i,j,k,l23);
		m13(i,j,k) = dlt2_sq * fsra * fsr[5] - dlt1_sq * Filter::TestFilter(i,j,k,l13);
	}}}

	/***** solve Lij *****/

	// Ui, cell-centered, stored in Lij
	#pragma omp single
	{
		u.Ugrid2CellCenter(uc);
		v.Vgrid2CellCenter(vc);
		w.Wgrid2CellCenter(wc);
	}

	// F(Ui) and Ui*Uj
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		// F(Ui), cell-centered, stored in Lii
		uf(i,j,k) = Filter::TestFilter(i,j,k,uc);
		vf(i,j,k) = Filter::TestFilter(i,j,k,vc);
		wf(i,j,k) = Filter::TestFilter(i,j,k,wc);
		// Ui*Uj
		s11(i,j,k) = uc(i,j,k) * uc(i,j,k);
		s22(i,j,k) = vc(i,j,k) * vc(i,j,k);
		s33(i,j,k) = wc(i,j,k) * wc(i,j,k);
		s12(i,j,k) = uc(i,j,k) * vc(i,j,k);
		s23(i,j,k) = vc(i,j,k) * wc(i,j,k);
		s13(i,j,k) = uc(i,j,k) * wc(i,j,k);
	}}}

	// Lij
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		// Lij = F(Ui) * F(Uj) - F(Ui*Uj)
		l12(i,j,k) = uf(i,j,k) * vf(i,j,k) - Filter::TestFilter(i,j,k,s12);
		l23(i,j,k) = vf(i,j,k) * wf(i,j,k) - Filter::TestFilter(i,j,k,s23);
		l13(i,j,k) = uf(i,j,k) * wf(i,j,k) - Filter::TestFilter(i,j,k,s13);
		// Lii = F(Ui)^2 - F(Ui^2), note: this must be after Lij
		l11(i,j,k) = uf(i,j,k) * uf(i,j,k) - Filter::TestFilter(i,j,k,s11);
		l22(i,j,k) = vf(i,j,k) * vf(i,j,k) - Filter::TestFilter(i,j,k,s22);
		l33(i,j,k) = wf(i,j,k) * wf(i,j,k) - Filter::TestFilter(i,j,k,s33);
		// Lij^d
		double iso = 1./3 * (l11(i,j,k) + l22(i,j,k) + l33(i,j,k));
		l11(i,j,k) -= iso;
		l22(i,j,k) -= iso;
		l33(i,j,k) -= iso;
	}}}

	/***** calculate eddy viscosity *****/

	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int i=1; i<ms.Nx; i++) {

		double lm = 0;
		double mm = 0;

		for (int k=1; k<ms.Nz; k++) {

			lm += l11(i,j,k) * m11(i,j,k)
			    + l22(i,j,k) * m22(i,j,k)
			    + l33(i,j,k) * m33(i,j,k)
			    + l12(i,j,k) * m12(i,j,k) * 2.
			    + l23(i,j,k) * m23(i,j,k) * 2.
			    + l13(i,j,k) * m13(i,j,k) * 2.;

			mm += m11(i,j,k) * m11(i,j,k)
			    + m22(i,j,k) * m22(i,j,k)
			    + m33(i,j,k) * m33(i,j,k)
			    + m12(i,j,k) * m12(i,j,k) * 2.
			    + m23(i,j,k) * m23(i,j,k) * 2.
			    + m13(i,j,k) * m13(i,j,k) * 2.;
		}

		for (int k=1; k<ms.Nz; k++)
			nut(i,j,k) *= fmin(fmax(.5*lm/mm, 0), .5) * pow(ms.dx(i)*ms.dy(j)*ms.dz(k), 2./3.); // |Sij| has been assigned to nut before
	}}
}

	// set boundary eddy viscosity
	Bcond::SetBoundaryY(nut, 1); // homogeneous Neumann
	Bcond::SetBoundaryX(nut, 3); // periodic
	Bcond::SetBoundaryZ(nut, 3); // periodic
}


void SGS::DynamicVreman(Scla &nut, const Vctr &vel, double Re)
{
	const Mesh &ms = nut.ms;
	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	// Scla &g11 = m11, &g12 = m22, &g13 = m33;
	// Scla &g21 = m12, &g22 = m23, &g23 = m13;
	// Scla &g31 = l11, &g32 = l22, &g33 = l33;
	// Scla &pss = l12, &gg  = l23;

	Scla s11(ms), s22(ms), s33(ms);
	Scla s12(ms), s23(ms), s13(ms);
	Scla g11(ms), g12(ms), g13(ms);
	Scla g21(ms), g22(ms), g23(ms);
	Scla g31(ms), g32(ms), g33(ms);
	Scla pss(ms), gg(ms);

	double sum1 = 0;
	double sum2 = 0;

#pragma omp parallel
{
	// Sij
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		const double *sr = vel.Strainrate(i,j,k);

		s11(i,j,k) = sr[0]; s22(i,j,k) = sr[1]; s33(i,j,k) = sr[2];
		s12(i,j,k) = sr[3]; s23(i,j,k) = sr[4]; s13(i,j,k) = sr[5];
	}}}

	// Gij, G:G, PI^g * S:S
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		// Gij = dUj/dXi
		const double *gr = vel.Gradient(i,j,k);

		g11(i,j,k) = gr[0]; g12(i,j,k) = gr[1]; g13(i,j,k) = gr[2];
		g21(i,j,k) = gr[3]; g22(i,j,k) = gr[4]; g23(i,j,k) = gr[5];
		g31(i,j,k) = gr[6]; g32(i,j,k) = gr[7]; g33(i,j,k) = gr[8];

		// G:G = (dUj/dXi)(dUj/dXi)
		double gr2 = 0;
		for (int n=0; n<9; n++) gr2 += gr[n] * gr[n];
		gg(i,j,k) = gr2;

		// PI^g = sqrt( B / G:G )
		double dx2 = ms.dx(i) * ms.dx(i);
		double dy2 = ms.dy(j) * ms.dy(j);
		double dz2 = ms.dz(k) * ms.dz(k);
		double beta[6];

		beta[0] = dx2 * gr[0]*gr[0] + dy2 * gr[3]*gr[3] + dz2 * gr[6]*gr[6];
		beta[1] = dx2 * gr[1]*gr[1] + dy2 * gr[4]*gr[4] + dz2 * gr[7]*gr[7];
		beta[2] = dx2 * gr[2]*gr[2] + dy2 * gr[5]*gr[5] + dz2 * gr[8]*gr[8];
		beta[3] = dx2 * gr[0]*gr[1] + dy2 * gr[3]*gr[4] + dz2 * gr[6]*gr[7];
		beta[4] = dx2 * gr[1]*gr[2] + dy2 * gr[4]*gr[5] + dz2 * gr[7]*gr[8];
		beta[5] = dx2 * gr[0]*gr[2] + dy2 * gr[3]*gr[5] + dz2 * gr[6]*gr[8];

		double bt2 = beta[0]*beta[1] + beta[1]*beta[2] + beta[0]*beta[2]
		           - beta[3]*beta[3] - beta[4]*beta[4] - beta[5]*beta[5];

		nut(i,j,k) = sqrt(bt2 / gr2);

		// PI^g * S:S =  PI^g * SijSij
		double sr2 = pow(s11(i,j,k), 2.)
		           + pow(s22(i,j,k), 2.)
		           + pow(s33(i,j,k), 2.)
		           + pow(s12(i,j,k), 2.) * 2.
		           + pow(s23(i,j,k), 2.) * 2.
		           + pow(s13(i,j,k), 2.) * 2.;

		pss(i,j,k) = nut(i,j,k) * sr2;
	}}}

	// C_nu
	#pragma omp barrier
	#pragma omp for reduction(+: sum1, sum2)
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double dx2 = ms.dx(i) * ms.dx(i);
		double dy2 = ms.dy(j) * ms.dy(j);
		double dz2 = ms.dz(k) * ms.dz(k);

		// F(Gij)
		double fgr[9];

		fgr[0] = Filter::TestFilter(i,j,k,g11);
		fgr[1] = Filter::TestFilter(i,j,k,g12);
		fgr[2] = Filter::TestFilter(i,j,k,g13);
		fgr[3] = Filter::TestFilter(i,j,k,g21);
		fgr[4] = Filter::TestFilter(i,j,k,g22);
		fgr[5] = Filter::TestFilter(i,j,k,g23);
		fgr[6] = Filter::TestFilter(i,j,k,g31);
		fgr[7] = Filter::TestFilter(i,j,k,g32);
		fgr[8] = Filter::TestFilter(i,j,k,g33);

		// FG:FG = F(dUj/dXi)F(dUj/dXi)
		double fgfg = 0;
		for (int n=0; n<9; n++) fgfg += fgr[n] * fgr[n];

		// FS:FS = F(Sij)F(Sij)
		double fsfs = pow(Filter::TestFilter(i,j,k,s11), 2.)
		            + pow(Filter::TestFilter(i,j,k,s22), 2.)
		            + pow(Filter::TestFilter(i,j,k,s33), 2.)
		            + pow(Filter::TestFilter(i,j,k,s12), 2.) * 2.
		            + pow(Filter::TestFilter(i,j,k,s23), 2.) * 2.
		            + pow(Filter::TestFilter(i,j,k,s13), 2.) * 2.;

		// PI^t
		double beta[6];
		beta[0] = 4.*dx2 * fgr[0]*fgr[0] + dy2 * fgr[3]*fgr[3] + 4.*dz2 * fgr[6]*fgr[6];
		beta[1] = 4.*dx2 * fgr[1]*fgr[1] + dy2 * fgr[4]*fgr[4] + 4.*dz2 * fgr[7]*fgr[7];
		beta[2] = 4.*dx2 * fgr[2]*fgr[2] + dy2 * fgr[5]*fgr[5] + 4.*dz2 * fgr[8]*fgr[8];
		beta[3] = 4.*dx2 * fgr[0]*fgr[1] + dy2 * fgr[3]*fgr[4] + 4.*dz2 * fgr[6]*fgr[7];
		beta[4] = 4.*dx2 * fgr[1]*fgr[2] + dy2 * fgr[4]*fgr[5] + 4.*dz2 * fgr[7]*fgr[8];
		beta[5] = 4.*dx2 * fgr[0]*fgr[2] + dy2 * fgr[3]*fgr[5] + 4.*dz2 * fgr[6]*fgr[8];

		double bt2 = beta[0]*beta[1] + beta[1]*beta[2] + beta[0]*beta[2]
		           - beta[3]*beta[3] - beta[4]*beta[4] - beta[5]*beta[5];

		double weight = ms.dx(i) * ms.dy(j) * ms.dz(k);

		sum1 += weight * (Filter::TestFilter(i,j,k,gg) - fgfg);
		sum2 += weight * (Filter::TestFilter(i,j,k,pss) - sqrt(bt2 / fgfg) * fsfs);
	}}}
}

	nut *= fmin(fmax(-.5/Re*sum1/sum2, 0), .5); // PI^g has been assigned to nut

	// set boundary eddy viscosity
	Bcond::SetBoundaryY(nut, 1); // homogeneous Neumann
	Bcond::SetBoundaryX(nut, 3); // periodic
	Bcond::SetBoundaryZ(nut, 3); // periodic
}





#define DIRTY_TRICK_SGS_


void SGS::SubGridStress(Vctr &shear, Vctr &normal, const Vctr &veldns, double rsclx, double rsclu)
// calculate sub-grid-scale stress by direct filter of resolved velocity
// normal & shear stress both at cell-centers, all up to virtual boundary
{
	const Mesh &ms = shear.ms; // refer to SPARSE mesh on which sgs stress is deployed
	const Mesh &ms0 = veldns.ms; // refer to RESOLVED mesh on which velocity is deployed

	const Scla &u = veldns[1];
	const Scla &v = veldns[2];
	const Scla &w = veldns[3];

	Scla &tau11 = normal[1], &tau12 = shear[1];
	Scla &tau22 = normal[2], &tau23 = shear[2];
	Scla &tau33 = normal[3], &tau13 = shear[3];

	Scla uu(ms0), uv(ms0), &uc = uu; u.Ugrid2CellCenter(uc); // cell-center-interpolation are for cross terms
	Scla vv(ms0), vw(ms0), &vc = vv; v.Vgrid2CellCenter(vc); // virtual boundary use linear extrapolation
	Scla ww(ms0), uw(ms0), &wc = ww; w.Wgrid2CellCenter(wc); // with periodicity ignored

	// calculate cross terms at cell-centers
	(uv.Set(uc)) *= vc;
	(vw.Set(vc)) *= wc; vv *= vc;
	(uw.Set(uc)) *= wc; uu *= uc; ww *= wc;

	// solve stress at cell-centers
#ifndef DIRTY_TRICK_SGS_
	#pragma omp parallel for
	for (int j=0; j<=ms.Ny; j++) {
#else
	for (int j=1; j<ms.Ny; j+=ms.Ny-2) {
	#pragma omp parallel for
#endif
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		double x = ms.xc(i) * rsclx, dx = ms.dx(i) * rsclx;
		double z = ms.zc(k) * rsclx, dz = ms.dz(k) * rsclx;
		double y = Filter::WallRscl(ms.yc(j), rsclx), dy = 0;

		tau11(i,j,k) = pow(Filter::FilterNodeU(x,y,z,dx,dy,dz,u), 2.) - Filter::FilterNodeA(x,y,z,dx,dy,dz,uu);
		tau22(i,j,k) = pow(Filter::FilterNodeV(x,y,z,dx,dy,dz,v), 2.) - Filter::FilterNodeA(x,y,z,dx,dy,dz,vv);
		tau33(i,j,k) = pow(Filter::FilterNodeW(x,y,z,dx,dy,dz,w), 2.) - Filter::FilterNodeA(x,y,z,dx,dy,dz,ww);

		tau12(i,j,k) =
			Filter::FilterNodeU(x,y,z,dx,dy,dz,u) *
			Filter::FilterNodeV(x,y,z,dx,dy,dz,v) -
			Filter::FilterNodeA(x,y,z,dx,dy,dz,uv);

		tau23(i,j,k) =
			Filter::FilterNodeV(x,y,z,dx,dy,dz,v) *
			Filter::FilterNodeW(x,y,z,dx,dy,dz,w) -
			Filter::FilterNodeA(x,y,z,dx,dy,dz,vw);

		tau13(i,j,k) =
			Filter::FilterNodeU(x,y,z,dx,dy,dz,u) *
			Filter::FilterNodeW(x,y,z,dx,dy,dz,w) -
			Filter::FilterNodeA(x,y,z,dx,dy,dz,uw);

		// handle HALFMFU
		if (ms.Ny < 2/ms0.Ly*j) tau12(i,j,k) *= -1;
		if (ms.Ny < 2/ms0.Ly*j) tau23(i,j,k) *= -1;

		// rescale to match viscous scale
		tau11(i,j,k) *= rsclu*rsclu; tau12(i,j,k) *= rsclu*rsclu;
		tau22(i,j,k) *= rsclu*rsclu; tau23(i,j,k) *= rsclu*rsclu;
		tau33(i,j,k) *= rsclu*rsclu; tau13(i,j,k) *= rsclu*rsclu;
	}}}
}


void SGS::SubGridShearStress(Vctr &shear, const Vctr &veldns, double rsclx, double rsclu)
// calculate sgs shear stress on edges up to virtual boundary
// by direct filter of resolved velocity
{
	const Mesh &ms = shear.ms; // refer to SPARSE mesh on which sgs stress is deployed
	const Mesh &ms0 = veldns.ms; // refer to RESOLVED mesh on which velocity is deployed

	const Scla &u = veldns[1];
	const Scla &v = veldns[2];
	const Scla &w = veldns[3];

	Scla uv(ms0), &tau12 = shear[1], &uc = uv; u.Ugrid2CellCenter(uc); // cell-center-interpolation are for cross terms
	Scla vw(ms0), &tau23 = shear[2], &vc = vw; v.Vgrid2CellCenter(vc); // virtual boundary use linear extrapolation
	Scla uw(ms0), &tau13 = shear[3], &wc = uw; w.Wgrid2CellCenter(wc); // with periodicity ignored

	// calculate cross terms at cell-centers
	Scla temp(uc);
	uv *= vc;
	vw *= wc;
	uw *= temp;

	// solve shear stress on edges
	#pragma omp parallel for collapse(3)
#ifndef DIRTY_TRICK_SGS_
	for (int j=0; j<=ms.Ny; j++) {
#else
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
#endif
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		double x = ms.x (i) * rsclx, dx = ms.hx(i) * rsclx;
		double z = ms.zc(k) * rsclx, dz = ms.dz(k) * rsclx;
		double y = Filter::WallRscl(ms.y(j), rsclx), dy = 0;

		if (i>0 && j>0) tau12(i,j,k) = (
			Filter::FilterNodeU(x,y,z,dx,dy,dz,u) *
			Filter::FilterNodeV(x,y,z,dx,dy,dz,v) -
			Filter::FilterNodeA(x,y,z,dx,dy,dz,uv) ) * pow(rsclu, 2.);

		x = ms.xc(i) * rsclx; dx = ms.dx(i) * rsclx;
		z = ms.z (k) * rsclx; dz = ms.hz(k) * rsclx;
		y = Filter::WallRscl(ms.y(j), rsclx); dy = 0;

		if (j>0 && k>0) tau23(i,j,k) = (
			Filter::FilterNodeV(x,y,z,dx,dy,dz,v) *
			Filter::FilterNodeW(x,y,z,dx,dy,dz,w) -
			Filter::FilterNodeA(x,y,z,dx,dy,dz,vw) ) * pow(rsclu, 2.);

#ifndef DIRTY_TRICK_SGS_
		x = ms.x (i) * rsclx; dx = ms.hx(i) * rsclx;
		z = ms.z (k) * rsclx; dz = ms.hz(k) * rsclx;
		y = Filter::WallRscl(ms.yc(j), rsclx); dy = 0;

		if (i>0 && k>0) tau13(i,j,k) = (
			Filter::FilterNodeU(x,y,z,dx,dy,dz,u) *
			Filter::FilterNodeW(x,y,z,dx,dy,dz,w) -
			Filter::FilterNodeA(x,y,z,dx,dy,dz,uw) ) * pow(rsclu, 2.);
#endif

		// handle HALFMFU
		if (1+ms.Ny < 2/ms0.Ly*j) tau12(i,j,k) *= -1;
		if (1+ms.Ny < 2/ms0.Ly*j) tau23(i,j,k) *= -1;
	}}}
}


void SGS::SubGridNormalStress(Vctr &normal, const Vctr &veldns, double rsclx, double rsclu)
// calculate sgs normal stress on cell-centers up to virtual boundary
// by direct filter of resolved velocity
{
	const Mesh &ms = normal.ms; // refer to SPARSE mesh on which sgs stress is deployed
	const Mesh &ms0 = veldns.ms; // refer to RESOLVED mesh on which velocity is deployed

	const Scla &u = veldns[1];
	const Scla &v = veldns[2];
	const Scla &w = veldns[3];

	Scla uu(ms0), &tau11 = normal[1]; (uu.Set(u)) *= u;
	Scla vv(ms0), &tau22 = normal[2]; (vv.Set(v)) *= v;
	Scla ww(ms0), &tau33 = normal[3]; (ww.Set(w)) *= w;

	// solve shear stress on edges
#ifndef DIRTY_TRICK_SGS_
	#pragma omp parallel for
	for (int j=0; j<=ms.Ny; j++) {
#else
	for (int j=1; j<ms.Ny; j+=ms.Ny-2) {
	#pragma omp parallel for
#endif
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		double x = ms.xc(i) * rsclx, dx = ms.dx(i) * rsclx;
		double z = ms.zc(k) * rsclx, dz = ms.dz(k) * rsclx;
		double y = Filter::WallRscl(ms.yc(j), rsclx), dy = 0;

		tau11(i,j,k) = (
			pow(Filter::FilterNodeU(x,y,z,dx,dy,dz,u), 2.) -
			Filter::FilterNodeU(x,y,z,dx,dy,dz,uu) ) * pow(rsclu, 2.);

		tau22(i,j,k) = (
			pow(Filter::FilterNodeV(x,y,z,dx,dy,dz,v), 2.) -
			Filter::FilterNodeV(x,y,z,dx,dy,dz,vv) ) * pow(rsclu, 2.);

		tau33(i,j,k) = (
			pow(Filter::FilterNodeW(x,y,z,dx,dy,dz,w), 2.) -
			Filter::FilterNodeW(x,y,z,dx,dy,dz,ww) ) * pow(rsclu, 2.);
	}}}
}






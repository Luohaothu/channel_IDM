#include "DA.h"
#include "Matrix.h"
#include "Interp.h"
#include "Bcond.h"

using namespace std;


DA::DA(const Mesh &ms):
fb_  (ms),
fldh_(ms),
vele_(ms),
msk_ (ms)
{
	erro_ = 0;
	iter_ = 0;
	fb_   = 0;
	fldh_ = 0;
	vele_ = 0;
	msk_  = 0;
	SetMsk(msk_);
}


bool DA::GetExp(double time, const Vctr &vele)
// wrap operations feeding reference data into assimilator
// may extend to file IO in the future
{
	double interval = 5e-2;

	if (fmod(time+INFTSM/2., interval) > INFTSM) {
		return false; // no experiment data provided
	}
	else {
		Interp::InterpBulkU(vele_[1], vele[1]);
		Interp::InterpBulkV(vele_[2], vele[2]);
		Interp::InterpBulkW(vele_[3], vele[3]);
	}
	return true;
}

void DA::Reset()
{
	fb_ = iter_ = 0;
}

bool DA::IfIter(double en, const Vctr &vel)
{
	erro_ = urhs(fldh_, vel, vele_, msk_);

	if (en < 1) return erro_ > en;
	else        return iter_ < en - INFTSM;
}

const Vctr& DA::GetAsmForce(const Vctr &vel, const Flow &vis, double dt, double alpha)
{
	CalcAdjoint(fldh_, vel, vis, dt);
	CalcForce(fb_, fldh_.SeeVec(), alpha);
	iter_++;
	return fb_;
}

void DA::WriteLog(double time)
{
	FILE *fp = fopen("DALOG.dat", "a");
	fprintf(fp, "DALOG Time: %f\tIter: %d\tError: %f\n", time, iter_, erro_);
	fclose(fp);
}

void DA::WriteForce(const char *path, int tstep) const
{
	char str[32];
	sprintf(str, "FDAX%08i", tstep); fb_[1].FileIO(path, str, 'w');
	sprintf(str, "FDAY%08i", tstep); fb_[2].FileIO(path, str, 'w');
	sprintf(str, "FDAZ%08i", tstep); fb_[3].FileIO(path, str, 'w');
}


// ***** non-APIs ***** //

void DA::SetMsk(Vctr &msk)
{
	const Mesh &ms = msk.ms;

	int Nx = ms.Nx;
	int Ny = ms.Ny;
	int Nz = ms.Nz;

	// number of layers to be assimilated above boundary (which is never assimilated)
	int nlayers = Ny/2 + 1;

	for (int j=0; j<=Ny; j++) {
	for (int k=0; k<=Nz; k++) {
	for (int i=0; i<=Nx; i++) {

		bool in_domain = (0<i && i<Nx) && (0<j && j<Ny) && (0<k && k<Nz);

		msk[1](i,j,k) = in_domain && (i>1) && (Ny-nlayers <= j || j <= nlayers);
		msk[2](i,j,k) = in_domain && (j>1) && (Ny-nlayers <= j || j <= nlayers+1);
		msk[3](i,j,k) = in_domain && (k>1) && (Ny-nlayers <= j || j <= nlayers);
	}}}
}

void DA::CalcAdjoint(Flow &fldh, const Vctr &vel, const Flow &vis, double dt)
{
	// urhs was called in IfIter() with U^n+1
	// now U has been rolled back to step n

	Vctr &velh = fldh.GetVec();
	Scla &dp   = fldh.GetScl();

	#pragma omp parallel
	{
		getuh1(velh, vel, vis, dt);
		getuh2(velh, vel, vis, dt);
		getuh3(velh, vel, vis, dt);
	}

	rhsdp(dp, velh, dt);
	dp.fftxz();
	#pragma omp parallel
	getfdp(dp, 0.);
	dp.ifftxz();

	update(velh, dp, dt);

	// no-slip BC in Y-direction applied automatically
	Bcond::SetBoundaryX(velh);
	Bcond::SetBoundaryZ(velh);
}

void DA::CalcForce(Vctr &fb, const Vctr &velh, double alpha)
{
	const Mesh &ms = fb.ms;
	
	double lambda = 0;

	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		lambda = fmax(lambda, velh.Module(i,j,k));
	}}}

	lambda = alpha / lambda;

	// steepest descent search of F
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		fb[1](i,j,k) += lambda * velh[1](i,j,k);
		fb[2](i,j,k) += lambda * velh[2](i,j,k);
		fb[3](i,j,k) += lambda * velh[3](i,j,k);
	}}}
}


// ***** adjoint velocity computation ***** //

double DA::urhs(Flow &fldh, const Vctr &vel, const Vctr &vele, const Vctr &msk)
// compute the RHS (mbc=cbc=0 for homogeneous Dirichlet BC) of momentum equation
{
	const Mesh &ms = fldh.ms;

	Vctr &velh = fldh.GetVec();
	Scla &temp = fldh.GetScl();

	Scla &uh = velh[1];
	Scla &vh = velh[2];
	Scla &wh = velh[3];

	const Scla &u = vel[1], &ue = vele[1], &msk1 = msk[1];
	const Scla &v = vel[2], &ve = vele[2], &msk2 = msk[2];
	const Scla &w = vel[3], &we = vele[3], &msk3 = msk[3];

	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		int id = ms.idx(i,j,k);
		
		if (i>1) uh[id] = -2 * msk1[id] * (u[id] - ue[id]);
		if (j>1) vh[id] = -2 * msk2[id] * (v[id] - ve[id]);
		if (k>1) wh[id] = -2 * msk3[id] * (w[id] - we[id]);

		temp(i,j,k) = msk.Module(i,j,k);
	}}}

	// return relative error, urhs need to be rescaled by 2
	double weight = .5 / temp.MeanA();

	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		temp(i,j,k) = velh.Module(i,j,k) / vele.Module(i,j,k);
	}}}

	return weight * temp.MeanA();
}


void DA::getuh1(Vctr &velh, const Vctr &vel, const Flow &vis, double dt)
/* compute deltaU^**, result returned by uh (the RHS ruh should be pre stored in uh ) */
{
	double u1, u2, v1, v2, w1, w2;
	double vis1, vis2, vis3, vis4, vis5, vis6;

	const Mesh &ms = velh.ms;

	const Scla &u = vel[1], &nuc = vis.SeeScl();
	const Scla &v = vel[2], &nuy = vis.SeeVec(2);
	const Scla &w = vel[3], &nuz = vis.SeeVec(3);

	Scla &uh = velh[1];
	Scla &vh = velh[2];
	Scla &wh = velh[3];

	double *ap = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ac = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *am = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ar = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];

	Matrix matx(ms.Nx-1);
	Matrix maty(ms.Ny-1);
	Matrix matz(ms.Nz-1);
	
	// ( I + dt M_11^2 )
	#pragma omp barrier
	#pragma omp for
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		for (int j=1; j<ms.Ny; j++) {

			int id   = ms.idx(i,j,k);
			int im   = ms.idx(ms.ima(i),j,k);
			int jp   = ms.idx(i,ms.jpa(j),k);
			int imjp = ms.idx(ms.ima(i),ms.jpa(j),k);

			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
			double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
			double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

			bool ifb3 = (j==1);
			bool ifb4 = (j==ms.Ny-1);

			v1 = .5/hxc * (v[id]*dxm + v[im]  *dxc);
			v2 = .5/hxc * (v[jp]*dxm + v[imjp]*dxc);
			vis3 = nuz[id];
			vis4 = nuz[jp];

			ap[j] = (-.25/hyp * (v1+v2)                       - 1./dyc * vis4/hyp             ) * dt;
			ac[j] = (-.25/dyc * (v1+v2) * (dyp/hyp - dym/hyc) + 1./dyc *(vis4/hyp + vis3/hyc) ) * dt + 1;
			am[j] = ( .25/hyc * (v1+v2)                       - 1./dyc * vis3/hyc             ) * dt;

			ar[j] = dt * uh[id];
		}

		maty.tdma(&am[1], &ac[1], &ap[1], &ar[1]); // apj at j=Ny-1 and amj at j=1 are redundant in tdma, thus no need for explicit removal
		
		for (int j=1; j<ms.Ny; j++) uh(i,j,k) = ar[j];
	}}

	// ( I + dt M_11^1 )
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			int id = ms.idx(i,j,k);
			int im = ms.idx(ms.ima(i),j,k);
			int ip = ms.idx(ms.ipa(i),j,k);

			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);

			vis1 = nuc[im];
			vis2 = nuc[id];

			ap[i] = (-u[id]/hxc - 2./hxc * vis2/dxc             ) * dt;
			ac[i] = (             2./hxc *(vis2/dxc + vis1/dxm) ) * dt + 1;
			am[i] = ( u[id]/hxc - 2./hxc * vis1/dxm             ) * dt;

			ar[i] = uh[id];
		}

		matx.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int i=1; i<ms.Nx; i++) uh(i,j,k) = ar[i];
	}}

	// ( I + dt M_11^3 )
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int i=1; i<ms.Nx; i++) {
		for (int k=1; k<ms.Nz; k++) {

			int id   = ms.idx(i,j,k);
			int im   = ms.idx(ms.ima(i),j,k);
			int kp   = ms.idx(i,j,ms.kpa(k));
			int imkp = ms.idx(ms.ima(i),j,ms.kpa(k));

			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
			double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
			double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

			w1 = .5/hxc * (w[id]*dxm + w[im]  *dxc);
			w2 = .5/hxc * (w[kp]*dxm + w[imkp]*dxc);
			vis5 = nuy[id];
			vis6 = nuy[kp];

			ap[k] = (-.25/hzp * (w1+w2)                       - 1./dzc * vis6/hzp             ) * dt;
			ac[k] = (-.25/dzc * (w1+w2) * (dzp/hzp - dzm/hzc) + 1./dzc *(vis6/hzp + vis5/hzc) ) * dt + 1;
			am[k] = ( .25/hzc * (w1+w2)                       - 1./dzc * vis5/hzc             ) * dt;

			ar[k] = uh[id];
		}

		matz.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);
		
		for (int k=1; k<ms.Nz; k++) uh(i,j,k) = ar[k];
	}}

	delete[] ap;
	delete[] ac;
	delete[] am;
	delete[] ar;
}


void DA::getuh2(Vctr &velh, const Vctr &vel, const Flow &vis, double dt)
/* compute deltaU^**, result returned by vh (the RHS rvh should be pre stored in uh ) */
{
	double u1, u2, v1, v2, w1, w2, m21uh;
	double vis1, vis2, vis3, vis4, vis5, vis6;

	const Mesh &ms = velh.ms;

	const Scla &u = vel[1], &nux = vis.SeeVec(1);
	const Scla &v = vel[2], &nuc = vis.SeeScl();
	const Scla &w = vel[3], &nuz = vis.SeeVec(3);

	Scla &uh = velh[1];
	Scla &vh = velh[2];
	Scla &wh = velh[3];

	double *ap = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ac = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *am = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ar = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];

	Matrix matx(ms.Nx-1);
	Matrix maty(ms.Ny-2);
	Matrix matz(ms.Nz-1);

	// ( I + dt M_22^2 )
	#pragma omp barrier
	#pragma omp for
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		for (int j=2; j<ms.Ny; j++) {

			int id   =      ms.idx(i,j,k);
			int ipjm =      ms.idx(ms.ipa(i),ms.jma(j),k);
			int ip, jp, kp; ms.ipx(i,j,k,ip,jp,kp);
			int im, jm, km; ms.imx(i,j,k,im,jm,km);

			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
			double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
			double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

			vis3 = nuc[jm];
			vis4 = nuc[id];

			ap[j] = (-v[id]/hyc - 2./hyc * vis4/dyc             ) * dt;
			ac[j] = (             2./hyc *(vis4/dyc + vis3/dym) ) * dt + 1;
			am[j] = ( v[id]/hyc - 2./hyc * vis3/dym             ) * dt;

			// m21uh
			double u0 = .25/hyc * ((u[id]+u[ip]) * dym + (u[jm]+u[ipjm]) * dyc);
			vis1 = nuz[id];
			vis2 = nuz[ip];
			u2 = .5 * (uh[id] + uh[ip]  );
			u1 = .5 * (uh[jm] + uh[ipjm]);

			m21uh = -u0/hyc * (u2-u1)
			        -1./dxc/hyc * (vis2 * (uh[ip]-uh[ipjm]) - vis1 * (uh[id]-uh[jm]));

			ar[j] = dt * (vh[id] - m21uh);
		}

		maty.tdma(&am[2], &ac[2], &ap[2], &ar[2]); // apj at j=Ny-1 and amj at j=2 are redundant in tdma, thus no need for explicit removal
		
		for (int j=2; j<ms.Ny; j++) vh(i,j,k) = ar[j];
	}}

	// ( I + dt M_22^1 )
	#pragma omp barrier
	#pragma omp for
	for (int j=2; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			int id   = ms.idx(i,j,k);
			int ip   = ms.idx(ms.ipa(i),j,k);
			int jm   = ms.idx(i,ms.jma(j),k);
			int ipjm = ms.idx(ms.ipa(i),ms.jma(j),k);

			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
			double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
			double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

			u2 = .5/hyc * (u[ip]*dym + u[ipjm]*dyc);
			u1 = .5/hyc * (u[id]*dym + u[jm]  *dyc);
			vis1 = nuz[id];
			vis2 = nuz[ip];

			ap[i] = (-.25/hxp * (u1+u2)                       - 1./dxc * vis2/hxp             ) * dt;
			ac[i] = (-.25/dxc * (u1+u2) * (dxp/hxp - dxm/hxc) + 1./dxc *(vis2/hxp + vis1/hxc) ) * dt + 1;
			am[i] = ( .25/hxc * (u1+u2)                       - 1./dxc * vis1/hxc             ) * dt;

			ar[i] = vh[id];
		}

		matx.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);
		
		for (int i=1; i<ms.Nx; i++) vh(i,j,k) = ar[i];
	}}

	// ( I + dt M_22^3 )
	#pragma omp barrier
	#pragma omp for
	for (int j=2; j<ms.Ny; j++) {
	for (int i=0; i<ms.Nx; i++) {
		for (int k=0; k<ms.Nz; k++) {

			int id   = ms.idx(i,j,k);
			int jm   = ms.idx(i,ms.jma(j),k);
			int kp   = ms.idx(i,j,ms.kpa(k));
			int jmkp = ms.idx(i,ms.jma(j),ms.kpa(k));

			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
			double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
			double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

			w2 = .5/hyc * (w[kp]*dym + w[jmkp]*dyc);
			w1 = .5/hyc * (w[id]*dym + w[jm]  *dyc);
			vis5 = nux[id];
			vis6 = nux[kp];

			ap[k] = (-.25/hzp * (w1+w2)                       - 1./dzc * vis6/hzp             ) * dt;
			ac[k] = (-.25/dzc * (w1+w2) * (dzp/hzp - dzm/hzc) + 1./dzc *(vis6/hzp + vis5/hzc) ) * dt + 1;
			am[k] = ( .25/hzc * (w1+w2)                       - 1./dzc * vis5/hzc             ) * dt;

			ar[k] = vh[id];
		}

		matz.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int k=1; k<ms.Nz; k++) vh(i,j,k) = ar[k];
	}}
		
	delete[] ap;
	delete[] ac;
	delete[] am;
	delete[] ar;
}


void DA::getuh3(Vctr &velh, const Vctr &vel, const Flow &vis, double dt)
/* compute deltaW^*, deltaV^*, deltaU^*, and update U^*, V^*, W^* (the RHS rwh should be pre stored in wh ) */
{
	double u1, u2, v1, v2, w1, w2;
	double m31uh, m32vh, m23wh, m12vh, m13wh;
	double vis1, vis2, vis3, vis4, vis5, vis6;

	const Mesh &ms = velh.ms;

	Scla &uh = velh[1];
	Scla &vh = velh[2];
	Scla &wh = velh[3];

	const Scla &u = vel[1], &nux = vis.SeeVec(1);
	const Scla &v = vel[2], &nuy = vis.SeeVec(2);
	const Scla &w = vel[3], &nuz = vis.SeeVec(3);
	const Scla              &nuc = vis.SeeScl();

	double *ap = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ac = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *am = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ar = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];

	Matrix matx(ms.Nx-1);
	Matrix maty(ms.Ny-1);
	Matrix matz(ms.Nz-1);

	// ( I + dt M_33^2 )
	#pragma omp barrier
	#pragma omp for
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		for (int j=1; j<ms.Ny; j++) {

			int id =              ms.idx(i,j,k);
			int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
			int im, jm, km;       ms.imx(i,j,k,im,jm,km);
			int ipjm, jpkm, ipkm; ms.pmx(i,j,k,ipjm,jpkm,ipkm);
			int imjm, jmkm, imkm; ms.mmx(i,j,k,imjm,jmkm,imkm);

			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
			double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
			double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

			bool ifb3 = (j==1);
			bool ifb4 = (j==ms.Ny-1);

			v2 = .5/hzc * (v[jp]*dzm + v[jpkm]*dzc);
			v1 = .5/hzc * (v[id]*dzm + v[km]  *dzc);
			vis3 = nux[id];
			vis4 = nux[jp];

			ap[j] = (-.25/hyp * (v1+v2)                       - 1./dyc * vis4/hyp             ) * dt;
			ac[j] = (-.25/dyc * (v1+v2) * (dyp/hyp - dym/hyc) + 1./dyc *(vis4/hyp + vis3/hyc) ) * dt + 1;
			am[j] = ( .25/hyc * (v1+v2)                       - 1./dyc * vis3/hyc             ) * dt;

			// m31uh
			double u0 = .25/hzc * ((u[id]+u[ip]) * dzm + (u[km]+u[ipkm]) * dzc);
			vis1 = nuy[id];
			vis2 = nuy[ip];
			u2 = .5 * (uh[id] + uh[ip]  );
			u1 = .5 * (uh[km] + uh[ipkm]);

			m31uh = -u0/hzc * (u2-u1)
			        -.1/dxc/hzc * (vis2 * (uh[ip]-uh[ipkm]) - vis1 * (uh[id]-uh[km]));

			// m32vh
			double v0 = .25/hzc * ((v[id]+v[jp]) * dzm + (v[km]+v[jpkm]) * dzc);
			vis3 = nux[id] * (1-ifb3);
			vis4 = nux[jp] * (1-ifb4);
			v2 = .5 * (vh[id] * (1-ifb3) + vh[jp]   * (1-ifb4));
			v1 = .5 * (vh[km] * (1-ifb3) + vh[jpkm] * (1-ifb4));

			m32vh = -v0/hzc * (v2-v1)
			        -1./dyc/hzc * (vis4 * (vh[jp]-vh[jpkm]) - vis3 * (vh[id]-vh[km]));

			ar[j] = dt * (wh[id] - m31uh - m32vh);
		}

		maty.tdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int j=1; j<ms.Ny; j++) wh(i,j,k) = ar[j];
	}}

	// ( I + dt M_33^1 )
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			int id   = ms.idx(i,j,k);
			int ip   = ms.idx(ms.ipa(i),j,k);
			int km   = ms.idx(i,j,ms.kma(k));
			int ipkm = ms.idx(ms.ipa(i),j,ms.kma(k));

			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
			double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
			double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

			u2 = .5/hzc * (u[ip]*dzm + u[ipkm]*dzc);
			u1 = .5/hzc * (u[id]*dzm + u[km]  *dzc);
			vis1 = nuy[id];
			vis2 = nuy[ip];

			ap[i] = (-.25/hxp * (u1+u2)                       - 1./dxc * vis2/hxp             ) * dt;
			ac[i] = (-.25/dxc * (u1+u2) * (dxp/hxp - dxm/hxc) + 1./dxc *(vis2/hxp + vis1/hxc) ) * dt + 1;
			am[i] = ( .25/hxc * (u1+u2)                       - 1./dxc * vis1/hxc             ) * dt;

			ar[i] = wh[id];
		}

		matx.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int i=1; i<ms.Nx; i++) wh(i,j,k) = ar[i];
	}}

	// ( I + dt M_33^3 )
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int i=1; i<ms.Nx; i++) {
		for (int k=1; k<ms.Nz; k++) {

			int id = ms.idx(i,j,k);
			int kp = ms.idx(i,j,ms.kpa(k));
			int km = ms.idx(i,j,ms.kma(k));

			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);

			vis5 = nuc[km];
			vis6 = nuc[id];

			ap[k] = (-w[id]/hzc - 2./hzc * vis6/dzc             ) * dt;
			ac[k] = (             2./hzc *(vis6/dzc + vis5/dzm) ) * dt + 1;
			am[k] = ( w[id]/hzc - 2./hzc * vis5/dzm             ) * dt;

			ar[k] = wh[id];
		}

		matz.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int k=1; k<ms.Nz; k++) wh(i,j,k) = ar[k];
	}}

	delete[] ap;
	delete[] ac;
	delete[] am;
	delete[] ar;

	// update dvh
	#pragma omp barrier
	#pragma omp for
	for (int j=2; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		int id   = ms.idx(i,j,k);
		int jm   = ms.idx(i,ms.jma(j),k);
		int km   = ms.idx(i,j,ms.kma(k));
		int kp   = ms.idx(i,j,ms.kpa(k));
		int jmkp = ms.idx(i,ms.jma(j),ms.kpa(k));

		double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
		double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
		double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
		double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
		double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

		// m23wh
		double w0 = .25/hyc * ((w[id]+w[kp]) * dym + (w[jm]+w[jmkp]) * dyc);
		vis5 = nux[id];
		vis6 = nux[kp];
		w2 = .5 * (wh[id] + wh[kp]  );
		w1 = .5 * (wh[jm] + wh[jmkp]);

		m23wh = -w0/hyc * (w2-w1)
		        -1./dzc/hyc * (vis6 * (wh[kp]-wh[jmkp]) - vis5 * (wh[id]-wh[jm]));

		vh[id] -= dt * m23wh;
	}}}
	
	// update duh
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		int id =              ms.idx(i,j,k);
		int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
		int im, jm, km;       ms.imx(i,j,k,im,jm,km);
		int imjp, jmkp, imkp; ms.mpx(i,j,k,imjp,jmkp,imkp);

		double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
		double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
		double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
		double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
		double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

		bool ifb3 = (j==1);
		bool ifb4 = (j==ms.Ny-1);

		// m12vh
		double v0 = .25/hxc * ((v[id]+v[jp]) * dxm + (v[im]+v[imjp]) * dzc);
		vis3 = nuz[id] * (1-ifb3);
		vis4 = nuz[jp] * (1-ifb4);
		v2 = .5 * (vh[id] * (1-ifb3) + vh[jp]   * (1-ifb4));
		v1 = .5 * (vh[im] * (1-ifb3) + vh[imjp] * (1-ifb4));

		m12vh = -v0/hxc * (v2-v1)
		        -1./dyc/hxc * (vis4 * (vh[jp]-vh[imjp]) - vis3 * (vh[id]-vh[im]));

		// m13wh
		double w0 = .25/hxc * ((w[id]+w[kp]) * dxm + (w[im]+w[imkp]) * dxc);
		vis5 = nuy[id];
		vis6 = nuy[kp];
		w2 = .5 * (wh[id] + wh[kp]  );
		w1 = .5 * (wh[im] + wh[imkp]);

		m13wh = -w0/hxc * (w2-w1)
		        -1./dzc/hxc * (vis6 * (wh[kp]-wh[imkp]) - vis5 * (wh[id]-wh[im]));

		uh[id] -= dt * (m12vh + m13wh);
	}}}
}



/***** projector computation *****/

void DA::rhsdp(Scla &rdp, const Vctr &velh, double dt)
/* compute RHS of Poisson equation and store in dp */
{
	const Mesh &ms = rdp.ms;
	const Scla &uh = velh[1];
	const Scla &vh = velh[2];
	const Scla &wh = velh[3];

	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		int id =              ms.idx(i,j,k);
		int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
		double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);

		bool ifb3 = (j==1);
		bool ifb4 = (j==ms.Ny-1);

		// ( Du^h - cbc ) / dt
		rdp[id] = 1./dt * (
			1./dxc * ( uh[ip] - uh[id] )
		+	1./dyc * ( vh[jp] * (1-ifb4) // mbc, cbc = 0 for homogeneous BC
			         - vh[id] * (1-ifb3) )
		+	1./dzc * ( wh[kp] - wh[id] ) );
	}}}
}

void DA::getfdp(Scla &fdp, double refp)
/* compute FDP (in Fourier space), result returned by fdp (the RHS frdp should be pre stored in fdp ) */
{
	const Mesh &ms = fdp.ms;

	double *ap = new double[ms.Ny];
	double *ac = new double[ms.Ny];
	double *am = new double[ms.Ny];
	double *ar = new double[ms.Ny];
	double *ai = new double[ms.Ny];

	Matrix maty(ms.Ny-1);

	#pragma omp for
	for (int k=0; k<ms.Nz-1; k++) { // i, k begin from 0 in fourier space
	for (int i=0; i<ms.Nxc; i++) { // negative k_x need not solve
		for (int j=1; j<ms.Ny; j++) {

			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
			double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

			bool ifb3 = (j==1);
			bool ifb4 = (j==ms.Ny-1);

			ap[j] = 1./dyc *   (1-ifb4) / hyp;
			ac[j] =-1./dyc * ( (1-ifb4) / hyp
				             + (1-ifb3) / hyc ) - ms.kx2(i) - ms.kz2(k);
			am[j] = 1./dyc *   (1-ifb3) / hyc;

			ar[j] = fdp[ms.idfxz(2*i  ,j,k)];
			ai[j] = fdp[ms.idfxz(2*i+1,j,k)];
		}

		// set reference pressure P(kx=0,kz=0,j=1) to 0
		// fdp(kx=0,kz=0,j=1) = -Nxz * fp(kx=0,kz=0,j=1) ==> <dp(j=1)> = - <p(j=1)>
		if (k==0 && i==0) {
			ap[1] = 0;
			ac[1] = 1;
			am[1] = 0;
			ar[1] = - (ms.Nx-1)*(ms.Nz-1) * refp;
			ai[1] = 0;
		}

		maty.tdma(&am[1], &ac[1], &ap[1], &ar[1]);
		maty.tdma(&am[1], &ac[1], &ap[1], &ai[1]);

		for (int j=1; j<ms.Ny; j++) {

			fdp[ms.idfxz(2*i  ,j,k)] = ar[j];
			fdp[ms.idfxz(2*i+1,j,k)] = ai[j];
		}
	}}
	/* note:
		In Huang's program, imag part of FDP at wavenumbers
		(kx=0,kz=0), (kx=0,kz=Nz/2), (kx=Nx/2,kz=0), (kx=Nx/2,kz=Nz/2)
		are explicitly set to 0. This is not needed here, since
		zero imag FRDP lead directly to zero imag FDP	*/

	delete[] ap;
	delete[] ac;
	delete[] am;
	delete[] ar;
	delete[] ai;
}

void DA::update(Vctr &velh, const Scla &dp, double dt)
/* project from U^* to U^n+1 using DP */
{
	const Mesh &ms = velh.ms;
	Scla &uh = velh[1];
	Scla &vh = velh[2];
	Scla &wh = velh[3];

	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		int id =              ms.idx(i,j,k);
		int im, jm, km;       ms.imx(i,j,k,im,jm,km);
		double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);

		if (i>0) uh[id] -= dt/hxc * (dp[id] - dp[im]);
		if (j>1) vh[id] -= dt/hyc * (dp[id] - dp[jm]);
		if (k>0) wh[id] -= dt/hzc * (dp[id] - dp[km]);
	}}}
}




#include "Solver.h"
#include "Interp.h"
#include "Filter.h"
#include "IDM.h"
#include "SGS.h"
#include "DA.h"

using namespace std;


Solver::Solver(const Mesh &ms):
ms   (ms),
fld_ (ms),
fldh_(ms),
vis_ (ms),
fb_  (ms),
bc_  (ms),
sbc_ (ms),
mpg_ {0},
time_(0)
{}



// ***** initiation ***** //

void Solver::init(double energy)
{
	InitChannel(fld_, bc_, sbc_, energy);
}

void Solver::init(const Flow &fld)
{
	InitFrom(fld_, fld);
}

// ***** default cumputation ***** //

void Solver::evolve(double Re, double dt, int sgstyp)
{
	set_time(time_+dt);

	CalcVis(vis_, get_vel(), Re, sgstyp);
	
	Bcond::ChannelNoSlip(bc_, sbc_, ms);

	CalcFb(fb_, mpg_);

	IDM::calc(fld_, fldh_, vis_, fb_, bc_, sbc_, dt);

	CalcMpg(mpg_, get_vel(), get_velh(), dt);

	Bcond::SetBoundaryY(get_vel(), bc_, sbc_);
	Bcond::SetBoundaryX(get_vel());
	Bcond::SetBoundaryZ(get_vel());
}

// ***** test computation ***** //

void Solver::evolve(double Re, double dt, int sgstyp, Solver &solver0)
{
	set_time(time_+dt);

	CalcVis(vis_, get_vel(), Re, sgstyp);
	ModifyBoundaryVis(vis_, get_vel(), solver0.get_vel(), Re);

	Bcond::ChannelDirichlet(bc_, sbc_, ms, solver0.get_vel());

	CalcFb(fb_, mpg_);

	IDM::calc(fld_, fldh_, vis_, fb_, bc_, sbc_, dt);

	CalcMpg(mpg_, get_vel(), get_velh(), dt, solver0.get_mpg());

	Bcond::SetBoundaryY(get_vel(), bc_, sbc_);
	Bcond::SetBoundaryX(get_vel());
	Bcond::SetBoundaryZ(get_vel());

	// Assimilate(solver0.get_vel(), dt, 5, .1);
}



// ***** non API ***** //

void Solver::InitChannel(Flow &fld, Boundaries &bc, Boundaries &sbc, double energy)
// initiate flow field (U, V, W, P, including all boundaries) from laminar with random fluctions
{
	const Mesh &ms = fld.ms;

	Scla &u = fld.GetVec(1);
	Scla &v = fld.GetVec(2);
	Scla &w = fld.GetVec(3);

	fld.InitRand(energy);

	// impose parabolic profile from laminar flow
	for (int j=1; j<ms.Ny; j++) {
		double ytilde = ms.yc(j) / ms.Ly;
		u.AddLyr(6 * ytilde * (1. - ytilde), j);
	}

	// modify flow rate
	u /= u.MeanU(); // bulk mean U = 1.0 due to non-dimensionalization
	w -= w.MeanW(); // note: boundaries are modified through using bulk functions

	Bcond::ChannelNoSlip(bc, sbc, ms);
	Bcond::SetBoundaryX(fld.GetVec());
	Bcond::SetBoundaryZ(fld.GetVec());
	Bcond::SetBoundaryY(fld.GetVec(), bc, sbc);
}

void Solver::InitFrom(Flow &fld, const Flow &fld0)
// initiate flow field (U, V, W, P, including all boundaries) from interpolation of existing fields
{
	const Mesh &ms = fld.ms;
	Scla &u = fld.GetVec(1);
	Scla &v = fld.GetVec(2);
	Scla &w = fld.GetVec(3);
	Scla &p = fld.GetScl();

	Interp::InterpBulkU(u, fld0.SeeVec(1));
	Interp::InterpBulkV(v, fld0.SeeVec(2));
	Interp::InterpBulkW(w, fld0.SeeVec(3));
	Interp::InterpBulkA(p, fld0.SeeScl());

	for (int j=0; j<=ms.Ny; j++) {
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {
		u(0,j,k) = p(0,j,k) = 0;
		v(i,0,k) = p(i,0,k) = 0;
		w(i,j,0) = p(i,j,0) = 0;
	}}}
}


void Solver::CalcVis(Flow &vis, const Vctr &vel, double Re, int sgstyp)
{
	Scla &nuc = (vis.GetScl() = 0);

	// if (2 <= sgstyp && sgstyp <= 4) {

	// 	static SGS sgs(vis.ms);

	// 	switch (sgstyp) {
	// 	case 2: sgs.Smargorinsky (nuc, vel, Re, .18); break;
	// 	case 3: sgs.DynamicSmarg (nuc, vel         ); break;
	// 	case 4: sgs.DynamicVreman(nuc, vel, Re     ); break;
	// 	}
	// }

	switch (sgstyp) {
	case 2: SGS::Smargorinsky (nuc, vel, Re, .18); break;
	case 3: SGS::DynamicSmarg (nuc, vel         ); break;
	case 4: SGS::DynamicVreman(nuc, vel, Re     ); break;
	}

	nuc += 1./Re;

	vis.CellCenter2Edge();
}

void Solver::ModifyBoundaryVis(Flow &vis, const Vctr &vel, double tau12)
{
	const Mesh &ms = vis.ms;

	Scla &nuz = vis.GetVec(3);

	const Scla &u = vel[1];
	const Scla &v = vel[2];

	for (int j=1; j<=ms.Ny; j+= ms.Ny-1) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double s12 = .5 * (
			1./ms.hx(i) * (v(i,j,k) - v(ms.ima(i),j,k)) +
			1./ms.hy(j) * (u(i,j,k) - u(i,ms.jma(j),k)) );

		nuz(i,j,k) = fmax(.5*tau12/s12, 0);
	}}}
}

void Solver::ModifyBoundaryVis(Flow &vis, const Vctr &vel, const Vctr &vel0, double Re)
{
	const Mesh &ms = vis.ms;
	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	Scla &nux = vis.GetVec(1);
	Scla &nuz = vis.GetVec(3);

	Vctr shear(ms);
	Vctr &normal = shear; // no intend to solve normal stress, it will be overwritten by shear stress

	SGS::SubGridStress(shear, normal, vel0);

	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
	for (int k=1; k<=ms.Nz; k++) {
	for (int i=1; i<=ms.Nx; i++) {
		// boundary sgs shear stress filtered from resolved velocity field
		double tau12 = shear[1](i,j,k);
		double tau23 = shear[2](i,j,k);

		// boundary strain rate of the coarse velocity field
		double s12 = .5 * (
			1./ms.hx(i) * (v(i,j,k) - v(ms.ima(i),j,k)) +
			1./ms.hy(j) * (u(i,j,k) - u(i,ms.jma(j),k)) );

		double s23 = .5 * (
			1./ms.hy(k) * (w(i,j,k) - w(i,ms.jma(j),k)) +
			1./ms.hz(j) * (v(i,j,k) - v(i,j,ms.kma(k))) );

		if (k<ms.Nz) nuz(i,j,k) = fmax(1./Re + tau12 / s12, 0);
		if (i<ms.Nx) nux(i,j,k) = fmax(1./Re + tau23 / s23, 0);
	}}}
}


void Solver::CalcFb(Vctr &fb, const double mpg[3])
{
	fb[1] = - mpg[0];
	fb[2] = - mpg[1];
	fb[3] = - mpg[2];
}
void Solver::CalcFb(Vctr &fb, const double mpg[3], const Vctr &f)
{
	CalcFb(fb, mpg);
	fb[1] += f[1];
	fb[2] += f[2];
	fb[3] += f[3];
}


void Solver::CalcMpg(double mpg[3], Vctr &vel, Vctr &velh, double dt, const double mpgref[3])
// solve the increment of mean pressure gradient at n+1/2 step
// given reference MFR at n+1 step or MPG at n+1/2 step
{
	const Mesh &ms = vel.ms;
	Scla &u = vel[1], &uh = velh[1];
	Scla &w = vel[3], &wh = velh[3];

	// solve mpg increment with streamwise flowrate 2.0 and spanwise 0
	double dmpg1 = mpgref ? mpgref[0]-mpg[0] : (u.MeanU()-1)/dt;
	double dmpg3 = mpgref ? mpgref[2]-mpg[2] :  w.MeanW()   /dt;

	// update mpg
	mpg[0] += dmpg1;
	mpg[2] += dmpg3;

	// complement mpg increment that was not included in the velocity update step
	for (int j=1; j<ms.Ny; j++) {

		u.MnsLyr(dt*dmpg1, j);
		w.MnsLyr(dt*dmpg3, j);

		uh.MnsLyr(dmpg1, j);
		wh.MnsLyr(dmpg3, j);
	}
}


void Solver::RollBack(Flow &fld, const Flow &fldh, double dt)
{
	Scla &u = fld.GetVec(1);
	Scla &v = fld.GetVec(2);
	Scla &w = fld.GetVec(3);
	Scla &p = fld.GetScl();

	const Scla &uh = fldh.SeeVec(1);
	const Scla &vh = fldh.SeeVec(2);
	const Scla &wh = fldh.SeeVec(3);
	const Scla &dp = fldh.SeeScl();

	( (u /= dt) -= uh ) *= dt;
	( (v /= dt) -= vh ) *= dt;
	( (w /= dt) -= wh ) *= dt;
	( (p /= dt) -= dp ) *= dt;
}


void Solver::RemoveSpanMean(Vctr &vel)
{
	const Mesh &ms = vel.ms;
	Scla &u = vel[1];
	Scla &v = vel[2];

	double um, *usm = new double[ms.Nx+1];
	double vm, *vsm = new double[ms.Nx+1];

	for (int j=1; j<ms.Ny; j++) {
		um = 0;
		for (int i=1; i<ms.Nx; i++)
			um += (usm[i]=u.MeanAz(i,j)) / (ms.Nx-1);

		usm[ms.Nx] = u.MeanAz(ms.Nz,j);

		for (int k=0; k<=ms.Nz; k++)
		for (int i=1; i<=ms.Nx; i++)
			u(i,j,k) -= usm[i] - um;
	}

	for (int j=2; j<ms.Ny; j++) {
		vm = 0;
		for (int i=1; i<ms.Nx; i++)
			vm += (vsm[i]=v.MeanAz(i,j)) / (ms.Nx-1);

		vsm[0]     = v.MeanAz(0,j);
		vsm[ms.Nx] = v.MeanAz(ms.Nx,j);

		for (int k=0; k<=ms.Nz; k++)
		for (int i=0; i<=ms.Nx; i++)
			v(i,j,k) -= vsm[i] - vm;
	}

	delete[] usm;
	delete[] vsm;
}


/***** data assimilation *****/

void Solver::Assimilate(const Vctr &velexp, double dt, double en, double a)
{
	static DA assimilator(ms);

	const Vctr &vel = fld_.SeeVec();
	
	if (assimilator.GetExp(time_, velexp)) {

		assimilator.Reset();

		// // store the unassimilated flow field for recovery later
		// Flow FLD_temp(mesh), FLDH_temp(mesh);
		// FLD_temp.V[1] = FLD.V[1]; FLDH_temp.V[1] = FLDH.V[1];
		// FLD_temp.V[2] = FLD.V[2]; FLDH_temp.V[2] = FLDH.V[2];
		// FLD_temp.V[3] = FLD.V[3]; FLDH_temp.V[3] = FLDH.V[3];
		// FLD_temp.S    = FLD.S;    FLDH_temp.S    = FLDH.S;


		while (assimilator.IfIter(en, vel)) { // converged or reached max iterations yet
			
			RollBack(fld_, fldh_, dt);        // roll back the flow fields to the old time step
			CalcFb(fb_, mpg_,                 // MPG cannot be rolled back, but it should not matter
				assimilator.GetAsmForce(vel, vis_, dt, a));   // compute adjoint state and apply assimilating force
			IDM::calc(fld_, fldh_, vis_, fb_, bc_, sbc_, dt); // solve the time step again under the new force
		}

		assimilator.WriteLog(time_);
		// assimilator.writeForce("", round(time/dt));

		// // recover the flow field to the unassimilated state
		// FLD.V[1] = FLD_temp.V[1]; FLDH.V[1] = FLDH_temp.V[1];
		// FLD.V[2] = FLD_temp.V[2]; FLDH.V[2] = FLDH_temp.V[2];
		// FLD.V[3] = FLD_temp.V[3]; FLDH.V[3] = FLDH_temp.V[3];
		// FLD.S    = FLD_temp.S;    FLDH.S    = FLDH_temp.S;
	}
}


void Solver::debug_Output(const char path[]) const
{
	int Ny = ms.Ny;

	// char path[1024] = "debug/";
	{
		char names[4][32] = {"U", "V", "W", "P"};
		fld_.SeeVec(1).debug_AsciiOutput(path, names[0], 0, Ny+1);
		fld_.SeeVec(2).debug_AsciiOutput(path, names[1], 1, Ny+1);
		fld_.SeeVec(3).debug_AsciiOutput(path, names[2], 0, Ny+1);
		fld_.SeeScl( ).debug_AsciiOutput(path, names[3], 0, Ny+1);
	}
	{
		char names[4][32] = {"UH", "VH", "WH", "DP"};
		fldh_.SeeVec(1).debug_AsciiOutput(path, names[0], 0, Ny+1);
		fldh_.SeeVec(2).debug_AsciiOutput(path, names[1], 1, Ny+1);
		fldh_.SeeVec(3).debug_AsciiOutput(path, names[2], 0, Ny+1);
		fldh_.SeeScl( ).debug_AsciiOutput(path, names[3], 0, Ny+1);
	}
	{
		char names[4][32] = {"NUX", "NUY", "NUZ", "NUC"};
		vis_.SeeVec(1).debug_AsciiOutput(path, names[0], 1, Ny+1);
		vis_.SeeVec(2).debug_AsciiOutput(path, names[1], 0, Ny+1);
		vis_.SeeVec(3).debug_AsciiOutput(path, names[2], 1, Ny+1);
		vis_.SeeScl( ).debug_AsciiOutput(path, names[3], 0, Ny+1);
	}
	{
		char names[3][32] = {"FBX", "FBY", "FBZ"};
		fb_[1].debug_AsciiOutput(path, names[0], 0, Ny+1);
		fb_[2].debug_AsciiOutput(path, names[1], 1, Ny+1);
		fb_[3].debug_AsciiOutput(path, names[2], 0, Ny+1);
	}
}






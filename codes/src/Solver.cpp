#include "Solver.h"
#include "Interp.h"
#include "Filter.h"
#include "IDM.h"
#include "SGS.h"
#include "WM.h"
#include "DA.h"
#include "PIO.h"

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
	// InitChannel(fld_, bc_, sbc_, energy);
	InitTbl(fld_, bc_, sbc_, energy);
}

void Solver::init(const Flow &fld)
{
	InitFrom(fld_, fld);
}

// ***** default cumputation ***** //


// #define TIME_TEST_SLV_

#ifndef TIME_TEST_SLV_

void Solver::evolve(double Re, double dt, int sgstyp)
{
	set_time(time_+dt);



	CalcVis(vis_, get_vel(), Re, sgstyp);

	// WM::LogLawWallShear(vis_, get_vel(), Re);

	// Bcond::ChannelNoSlip(bc_, sbc_, ms);
	// Bcond::ChannelRobin(bc_, sbc_, get_vel(), vis_, 1.);
	// Bcond::ChannelHalf(bc_, sbc_, ms);
	Bcond::TblDevelop(bc_, sbc_, get_vel(), 1., dt);

	CalcFb(fb_, mpg_);

	IDM::calc(fld_, fldh_, vis_, fb_, bc_, sbc_, dt);

	Bcond::SetBoundaryX(get_vel(), bc_, sbc_);
	Bcond::SetBoundaryY(get_vel(), bc_, sbc_);
	Bcond::SetBoundaryZ(get_vel());

	// CalcMpg(mpg_, get_vel(), get_velh(), dt);

	if (sgstyp == 1) RemoveSpanMean(get_vel());
}

#endif

// ***** test computation ***** //

void Solver::evolve(double Re, double dt, int sgstyp, Solver &solver0, double Re0)
{
	set_time(time_+dt);

	double utau  = sqrt(fabs(mpg_[0]));
	double utau0 = sqrt(fabs(solver0.get_mpg()[0]));
	double rsclu = utau0 ? utau/utau0 : 1;
	double rsclx = rsclu * Re / Re0; // note: tiled MFU must satisfy periodic condition to a high presicion
	double Ret = utau * Re;

	// cout << utau << endl;
	// cout << utau0 << endl;
	// cout << rsclu << endl; // utau_LES / utau_MFU
	// cout << rsclx << endl; // Ret_LES / Ret_MFU
	// cout << Ret << endl;

	CalcVis(vis_, get_vel(), Re, sgstyp);

	WM::OffWallSubGridUniform(vis_, get_vel(), solver0.get_vel(), Re, Ret, rsclx, rsclu);

	Bcond::ChannelDirichlet(bc_, sbc_, ms, PIO::BoundaryPredict(get_vel(), solver0.get_vel(), Ret, rsclx, rsclu));
	// Bcond::ChannelDirichlet(bc_, sbc_, ms, solver0.get_vel(), rsclx, rsclu);

	CalcFb(fb_, mpg_);
	// fb_[1] *= 2; // accerate convergence

	IDM::calc(fld_, fldh_, vis_, fb_, bc_, sbc_, dt);

	Bcond::SetBoundaryY(get_vel(), bc_, sbc_);
	Bcond::SetBoundaryX(get_vel());
	Bcond::SetBoundaryZ(get_vel());

	// CalcMpg(mpg_, get_vel(), get_velh(), dt, solver0.get_mpg());

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
}

void Solver::InitTbl(Flow &fld, Boundaries &bc, Boundaries &sbc, double energy)
// initiate flow field (U, V, W, P, including all boundaries) from laminar with random fluctions
{
	const Mesh &ms = fld.ms;

	Scla &u = fld.GetVec(1);
	Scla &v = fld.GetVec(2);
	Scla &w = fld.GetVec(3);

	// impose initial profile
	for (int j=0; j<=ms.Ny; j++)
		u.AddLyr(fmin(10*ms.yc(j), 1.), j);

	Bcond::TblDevelop(bc, sbc, fld.SeeVec(), 1., 1);

	Bcond::SetBoundaryX(fld.GetVec(), bc, sbc);
	Bcond::SetBoundaryY(fld.GetVec(), bc, sbc);
	Bcond::SetBoundaryZ(fld.GetVec());
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
	Scla &nuc = vis.GetScl();

	if (sgstyp <= 1) {
		nuc.Set(1./Re);
	}
	else {
		static SGS sgs(nuc.ms);

		switch (sgstyp) {
		case 2: sgs.Smargorinsky (nuc, vel, Re, .18); break;
		case 3: sgs.DynamicSmarg (nuc, vel         ); break;
		case 4: sgs.DynamicVreman(nuc, vel, Re     ); break;
		}

		nuc += 1./Re;
	}

	vis.CellCenter2Edge();
}


void Solver::CalcFb(Vctr &fb, const double mpg[3])
{
	fb[1].Set(- mpg[0]);
	fb[2].Set(- mpg[1]);
	fb[3].Set(- mpg[2]);
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
	#pragma omp parallel for
	for (int j=1; j<ms.Ny; j++) {
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		if (i>0) u(i,j,k) -= dt * dmpg1;
		if (k>0) w(i,j,k) -= dt * dmpg3;

		if (i>0) uh(i,j,k) -= dmpg1;
		if (k>0) wh(i,j,k) -= dmpg3;
	}}}

	// note: since BC has been applied before, this function should maintain boundary relations
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

#pragma omp parallel
{
	double um, *usm = new double[ms.Nx+1];
	double vm, *vsm = new double[ms.Nx+1];

	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
		um = 0;
		for (int i=1; i<ms.Nx; i++)
			um += (usm[i]=u.MeanAz(i,j)) / (ms.Nx-1);

		usm[ms.Nx] = u.MeanAz(ms.Nz,j);

		for (int k=0; k<=ms.Nz; k++)
		for (int i=1; i<=ms.Nx; i++)
			u(i,j,k) -= usm[i] - um;
	}

	#pragma omp for
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




#ifdef TIME_TEST_SLV_

#include <sys/time.h>
int cnt = 0;
struct timeval time0, time1;
double t_vis=0, t_bc=0, t_fb=0, t_idm=0, t_mpg=0, t_setbc=0;

void Solver::evolve(double Re, double dt, int sgstyp)
{
	set_time(time_+dt);

	cout << ++ cnt << endl;

	gettimeofday(&time0, NULL);	CalcVis(vis_, get_vel(), Re, sgstyp);		gettimeofday(&time1, NULL);	t_vis += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	gettimeofday(&time0, NULL);	Bcond::ChannelNoSlip(bc_, sbc_, ms);		gettimeofday(&time1, NULL);	t_bc  += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	gettimeofday(&time0, NULL);	CalcFb(fb_, mpg_);							gettimeofday(&time1, NULL);	t_fb  += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	gettimeofday(&time0, NULL);	IDM::calc(fld_,fldh_,vis_,fb_,bc_,sbc_,dt);	gettimeofday(&time1, NULL);	t_idm += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	gettimeofday(&time0, NULL);	CalcMpg(mpg_, get_vel(), get_velh(), dt);	gettimeofday(&time1, NULL);	t_mpg += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	gettimeofday(&time0, NULL);	Bcond::SetBoundaryY(get_vel(), bc_, sbc_);
								Bcond::SetBoundaryX(get_vel());
								Bcond::SetBoundaryZ(get_vel());				gettimeofday(&time1, NULL);	t_setbc+=(time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	if (cnt == 10) {
		cout << endl << endl
			 << "total:\t"  << 1./cnt * (t_vis+t_bc+t_fb+t_idm+t_mpg+t_setbc) << endl
			 << "t_vis:\t"  << t_vis  / t_idm << endl
			 << "t_bc:\t"   << t_bc   / t_idm << endl
			 << "t_fb:\t"   << t_fb   / t_idm << endl
			 << "t_mpg:\t"  << t_mpg  / t_idm << endl
			 << "t_setbc:\t"<< t_setbc/ t_idm << endl << endl;
		exit(0);
	}
}

#endif





// void Solver::ModifyBoundaryVis(Flow &vis, const Vctr &vel, double tau12)
// {
// 	const Mesh &ms = vis.ms;

// 	Scla &nuz = vis.GetVec(3);

// 	const Scla &u = vel[1];
// 	const Scla &v = vel[2];

// 	for (int j=1; j<=ms.Ny; j+= ms.Ny-1) {
// 	for (int k=1; k<ms.Nz; k++) {
// 	for (int i=1; i<ms.Nx; i++) {

// 		double s12 = .5 * (
// 			1./ms.hx(i) * (v(i,j,k) - v(ms.ima(i),j,k)) +
// 			1./ms.hy(j) * (u(i,j,k) - u(i,ms.jma(j),k)) );

// 		nuz(i,j,k) = fmax(.5*tau12/s12, 0);
// 	}}}
// }

// void Solver::ModifyBoundaryVis(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re)
// {
// 	const Mesh &ms = vis.ms, &ms0 = veldns.ms;

// 	const Scla &u = vel[1], &u0 = veldns[1];
// 	const Scla &v = vel[2], &v0 = veldns[2];
// 	const Scla &w = vel[3], &w0 = veldns[3];

// 	// viscosity & sgs-stress on edges aimed for
// 	Vctr shearsgs(ms);
// 	Scla &nux = vis.GetVec(1);
// 	Scla &nuz = vis.GetVec(3);
// 	const Scla &tau23sgs = shearsgs[2];
// 	const Scla &tau12sgs = shearsgs[1];

// 	// sgs shear stress filtered from resolved velocity field
// 	SGS::SubGridShearStress(shearsgs, veldns);

// 	// for top & bottom real boundaries
// 	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

// 		double r12dns = 0; // DNS Reynolds stress
// 		double r12les = 0; // LES resolved Reynolds
// 		double v12sgs = 0; // mean sgs-stress filtered from DNS

// 		// the wall normal position on which stresses act
// 		const double y = ms.y(j);

// 		// ***** rescale the DNS-filtered sgs stress by Reynolds stress defect ***** //

// 		// calculate reference (DNS) Reynolds shear stress
// 		int j3u = Interp::BiSearch(y, ms0.yc(), 0, ms0.Ny), j4u = j3u+1;
// 		int j3v = Interp::BiSearch(y, ms0.y(),  1, ms0.Ny), j4v = j3v+1;

// 		double um3 = u0.meanxz(j3u), y3u = ms0.yc(j3u);
// 		double um4 = u0.meanxz(j4u), y4u = ms0.yc(j4u);
// 		double vm3 = v0.meanxz(j3v), y3v = ms0.y(j3v);
// 		double vm4 = v0.meanxz(j4v), y4v = ms0.y(j4v);

// 		double um = (um3 * (y4u-y) + um4 * (y-y3u)) / (y4u-y3u);
// 		double vm = (vm3 * (y4v-y) + vm4 * (y-y3v)) / (y4v-y3v);

// 		#pragma omp parallel for reduction(+: r12dns)
// 		for (int k=1; k<ms0.Nz; k++) { double z = ms0.zc(k);
// 		for (int i=1; i<ms0.Nx; i++) { double x = ms0.x(i);
// 			// calculate uv on z-edges and interpolate to the desired y position
// 			r12dns += (Interp::InterpNodeU(x,y,z,u0) - um)
// 			        * (Interp::InterpNodeV(x,y,z,v0) - vm)
// 			        * ((ms.Nx-1) * (ms.Nz-1)) / ((ms0.Nx-1) * (ms0.Nz-1)); // rescale to match the other two
// 		}}

// 		// calculate resolved Reynolds stress and mean sgs-stress
// 		um3 = u.meanxz(j3u = ms.jma(j)); y3u = ms.yc(j3u);
// 		um4 = u.meanxz(j4u = j);         y4u = ms.yc(j4u);
// 		vm  = v.meanxz(j);
// 		um  = (um3 * (y4u-y) + um4 * (y-y3u)) / (y4u-y3u);

// 		for (int k=1; k<ms.Nz; k++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hyc=ms.hy(j);
// 		for (int i=1; i<ms.Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxc=ms.hx(i);

// 			int id =        ms.idx(i,j,k);
// 			int im, jm, km; ms.imx(i,j,k,im,jm,km);

// 			v12sgs += tau12sgs[id];
// 			r12les += (.5/hyc * (u[id]*dym + u[jm]*dyc) - um)
// 			        * (.5/hxc * (v[id]*dxm + v[im]*dxc) - vm);
// 		}}

// 		const double rescale1 = fabs((r12dns - r12les) / v12sgs);
// 		const double rescale2 = fabs((y4u-y3u) / (1-fabs(y-1)) / log((1-fabs(y4u-1))/(1-fabs(y3u-1))));

// 		// // rescale physical viscosity by replacing low-order differencing with log law
// 		// double dyU2 = (um4-um3) / log((1-fabs(y4u-1))/(1-fabs(y3u-1))) / (1-fabs(y-1));
// 		// double dyU1 = (um4-um3) / (y4u-y3u);
// 		// const double rescale2 = fabs(dyU2 / dyU1);

// 		// ***** rescale viscosity to account for low-order differencing error ***** //

// 		// construct low- & high-order differenced U gradient
// 		// int j5u = j==1 ? j4u+1 : j3u-1;
// 		// double um5 = u.meanxz(j5u);
// 		// double y5u = ms.yc(j5u);

// 		// double dyU2 = (y4u+y5u-2*y) / (y4u-y3u) / (y3u-y5u) * um3
// 		//             + (y3u+y5u-2*y) / (y3u-y4u) / (y4u-y5u) * um4
// 		//             + (y3u+y4u-2*y) / (y3u-y5u) / (y5u-y4u) * um5;

// 		cout << rescale1 << " " << rescale2 << endl;

// 		// ***** modify physical & eddy viscosity accordingly ***** //

// 		for (int k=1; k<=ms.Nz; k++) {
// 		for (int i=1; i<=ms.Nx; i++) {
// 			// boundary strain rate of the coarse velocity field
// 			const double *sr = vel.ShearStrain(i,j,k);

// 			double s12 = sr[0], tau12 = tau12sgs(i,j,k) * rescale1;
// 			double s23 = sr[1], tau23 = tau23sgs(i,j,k) * rescale1;

// 			if (k<ms.Nz) nuz(i,j,k) = rescale2/Re + fmax(.5 * tau12 / s12, 0);
// 			if (i<ms.Nx) nux(i,j,k) = rescale2/Re + fmax(.5 * tau23 / s23, 0);
// 		}}
// 	}
// }






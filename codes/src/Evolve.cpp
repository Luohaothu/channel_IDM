#include "Solver.h"
#include "Interp.h"
#include "Filter.h"
#include "IDM.h"
#include "SGS.h"
#include "WM.h"
#include "OFW.h"
#include "Statis.h"

using namespace std;


void Solver::Manipulation()
{
	// CalcMpg(mpg, fldh.GetVec(), para.dt);
	CalcMpg(mpg, fldh.GetVec(), fb, para.dt);

	if (para.bftype == 1)
		RemoveSpanMean(fldh.GetVec());
}

// ***** choose cumputation mode ***** //

// #define _FULLCHANNEL
// #define _HALFCHANNEL
// #define _FULLCHANNEL_ALGEWM
// #define _FULLCHANNEL_SLIPWM
// #define _TBL
#define _EBL
// #define _SYNMFU
#define _TESTCHANNEL

// *********************************** //

#ifdef _FULLCHANNEL // full channel with noslip BC
void Solver::Evolve()
{
	step ++; time += para.dt;

	CalcVis(vis, fld.SeeVec(), para.Re, para.bftype);
	Bcond::ChannelNoSlip(bc, sbc, ms);
	CalcFb(fb, mpg);
	IDM::calc(fldh, fld, vis, fb, bc, sbc, para.dt);
	SetBoundaries(fldh.GetVec(), bc, sbc);
	Manipulation();

	fld += (fldh -= fld);
}
#endif

#ifdef _HALFCHANNEL
void Solver::Evolve() // half channel
{
	step ++; time += para.dt;

	CalcVis(vis, fld.SeeVec(), para.Re, para.bftype);
	Bcond::ChannelHalf(bc, sbc, ms);
	CalcFb(fb, mpg);
	IDM::calc(fldh, fld, vis, fb, bc, sbc, para.dt);
	SetBoundaries(fldh.GetVec(), bc, sbc);
	Manipulation();

	fld += (fldh -= fld);
}
#endif

#ifdef _FULLCHANNEL_ALGEWM // full channel with algebra wall model
void Solver::Evolve()
{
	step ++; time += para.dt;

	CalcVis(vis, fld.SeeVec(), para.Re, para.bftype);
	WM::LogLawWallShear(vis, fld.SeeVec(), para.Re);
	Bcond::ChannelNoSlip(bc, sbc, ms);
	CalcFb(fb, mpg);
	IDM::calc(fldh, fld, vis, fb, bc, sbc, para.dt);
	SetBoundaries(fldh.GetVec(), bc, sbc);
	Manipulation();

	fld += (fldh -= fld);
}
#endif

#ifdef _FULLCHANNEL_SLIPWM // full channel with slip wall model
void Solver::Evolve()
{
	step ++; time += para.dt;

	CalcVis(vis, fld.SeeVec(), para.Re, para.bftype);
	Bcond::ChannelRobin(bc, sbc, fld.SeeVec(), vis, 1.);
	CalcFb(fb, mpg);
	IDM::calc(fldh, fld, vis, fb, bc, sbc, para.dt);
	SetBoundaries(fldh.GetVec(), bc, sbc);
	Manipulation();

	fld += (fldh -= fld);
}
#endif

#ifdef _TBL // boundary layer with Ufree = 1
void Solver::Evolve()
{
	step ++; time += para.dt;

	CalcVis(vis, fld.SeeVec(), para.Re, para.bftype);
	Bcond::TblDevelop(bc, sbc, fld.SeeVec(), 1., para.dt);
	CalcFb(fb, mpg);
	IDM::calc(fldh, fld, vis, fb, bc, sbc, para.dt);
	SetBoundaries(fldh.GetVec(), bc, sbc);
	Manipulation();

	fld += (fldh -= fld);
}
#endif

#ifdef _EBL
void Solver::Evolve()
{
	step ++; time += para.dt;

	if (step==1) {mpg[0] = -1.; mpg[1] = mpg[2] = 0;}

	CalcVis(vis, fld.SeeVec(), para.Re, para.bftype);
	Bcond::ChannelHalf(bc, sbc, ms);
	CalcFb(fb, mpg, "Fx.txt");
	// Bcond::TblEquiv(bc, sbc, fld.SeeVec(), fb, vis, para.dt);
	IDM::calc(fldh, fld, vis, fb, bc, sbc, para.dt);
	SetBoundaries(fldh.GetVec(), bc, sbc);
	Manipulation();

	fld += (fldh -= fld);
}
#endif

#ifdef _TESTCHANNEL // off-wall channel
void Solver::Evolve(const Solver &slv)
{
	step ++; time += para.dt;

	// mpg[0] = -1.;

	double Ret   = InnerScale();
	double rsclx = Ret / slv.InnerScale(); // Ret_LES  / Ret_MFU
	double rsclu = rsclx * slv.para.Re / para.Re; // utau_LES / utau_MFU
	
	CalcVis(vis, fld.SeeVec(), para.Re, para.bftype);
	OFW::OffWallSubGridUniform(vis, fld.SeeVec(), slv.fld.SeeVec(), para.Re, Ret, rsclx, rsclu);
	OFW::OffWallVelo(bc, sbc, fldh.GetVec(), fld.SeeVec(), slv.fld.SeeVec(), Ret, rsclx, rsclu);
	CalcFb(fb, mpg);
	IDM::calc(fldh, fld, vis, fb, bc, sbc, para.dt);
	SetBoundaries(fldh.GetVec(), bc, sbc);
	Manipulation();

	fld += (fldh -= fld);
}
#endif

#ifdef _SYNMFU // off-wall channel with synthetic MFU
void Solver::Evolve()
{
	step ++; time += para.dt;
	
	CalcVis(vis, fld.SeeVec(), para.Re, para.bftype);
	WM::UniformReyShear(vis, fld.SeeVec(), -.9419, 0, para.Re);
	OFW::OffWallVelo(bc, sbc, fldh.GetVec(), fld.SeeVec(), time, InnerScale());
	CalcFb(fb, mpg);
	IDM::calc(fldh, fld, vis, fb, bc, sbc, para.dt);
	SetBoundaries(fldh.GetVec(), bc, sbc);
	Manipulation();

	fld += (fldh -= fld);
}
#endif




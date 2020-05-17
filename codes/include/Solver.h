#pragma once

#include "Basic.h"
#include "Bcond.h"


class Solver
{
public:
	const Mesh &ms;

	Solver(const Mesh &ms);

	// ***** initiation ***** //
	void init(double energy);
	void init(const Flow &fld);

	// ***** default cumputation ***** //
	void evolve(double Re, double dt, int sgstyp);

	// ***** test computation ***** //
	void evolve(double Re, double dt, int sgstyp, Solver &solver0);


	double get_time()   const { return time_; };
	double set_time(double t) { return time_ = t; };

	Flow& get_fld () { return fld_; };
	Flow& get_fldh() { return fldh_; };
	Flow& get_vis () { return vis_; };
	Vctr& get_fb  () { return fb_; };
	Vctr& get_vel () { return fld_.GetVec(); };
	Vctr& get_velh() { return fldh_.GetVec(); };
	double* get_mpg() { return mpg_; };

	void debug_Output(const char path[]) const;

private:
	Flow fld_;       // solution field: Velocity & Pressure
	Flow fldh_;      // time derivative of solution field, also serve as intermediate fields
	Flow vis_;       // viscous field: deployed on cell-centers and edges
	Vctr fb_;        // body force
	Boundaries bc_;  // boundary field specifying BCs at n+1 step
	Boundaries sbc_; // secondary boundary coefficients determined by the type of BC used
	double mpg_[3];  // mean pressure gradient acting as drving force
	double time_;

	// ***** initiation ***** //
	static void InitFrom(Flow &fld, const Flow &fld0);
	static void InitChannel(Flow &fld, Boundaries &bc, Boundaries &sbc, double energy);

	// ***** construct viscosity ***** //
	static void CalcVis(Flow &vis, const Vctr &vel, double Re, int sgstyp);
	static void ModifyBoundaryVis(Flow &vis, const Vctr &vel, double tau12);
	static void ModifyBoundaryVis(Flow &vis, const Vctr &vel, const Vctr &vel0, double Re);

	// ***** construct body forces ***** //
	static void CalcFb(Vctr &fb, const double mpg[3]);
	static void AddFb(Vctr &fb, const Vctr &f);

	// ***** adjusting mean pressure gradient ***** //
	static void CalcMpg(double mpg[3], Vctr &vel, Vctr &velh, double dt, const double mpgref[3] = NULL);
	
	// ***** data manipulation ***** //
	static void RollBack(Flow &fld, const Flow &fldh, double dt);
	static void RemoveSpanMean(Vctr &vel);
	
	void Assimilate(const Vctr &velexp, double dt, double en, double a);

	// ***** construct boundary conditions ***** //

	// Bcond::ChannelNoSlip(bc_, sbc_, ms);
	// Bcond::ChannelDirichlet(bc_, sbc_, ms, solver0.get_fld().SeeVec());
	// Bcond::ChannelRobin(bc_, sbc_, ms, solver0.get_fld().SeeVec());
	// Bcond::TblCycling(bc_, sbc_, ms, double ufree);

	// ***** apply boundary conditions ***** //

	// non-periodic
	// Bcond::SetBoundaryX(fld_.GetVec(), bc_, sbc_);
	// Bcond::SetBoundaryY(fld_.GetVec(), bc_, sbc_);
	// periodic
	// Bcond::SetBoundaryX(fld_.GetVec());
	// Bcond::SetBoundaryZ(fld_.GetVec());

};




		




#pragma once

#include "Basic.h"
#include "Bcond.h"
#include "Para.h"


class Solver
{
public:
	Solver(const char* path, bool ifinit=true);

	void ContinueCase();
	void Manipulation();

	void Evolve(); // default computation
	void Evolve(const Solver &slv); // test computation

	void Output();
	void ShowInfo() const;
	void debug_Output(const char path[]) const;

	Para para;
	Geometry geom;
	const Mesh ms;

	Boundaries bc;  // boundary field specifying BCs at n+1 step
	Boundaries sbc; // secondary boundary coefficients determined by the type of BC used
	Flow fld;       // solution field: Velocity & Pressure
	Flow fldh;      // time derivative of solution field, also serve as intermediate fields
	Flow vis;       // viscous field: deployed on cell-centers and edges
	Vctr fb;        // body force
	std::vector<double> mpg; // mean pressure gradient acting as drving force
	
	double time;
	int step;

private:
	// ***** initiation ***** //
	void InitFieldFrom(const Flow &fld0);
	void InitFieldChan(Boundaries &bc, Boundaries &sbc, double energy);
	void InitFieldBlyr(Boundaries &bc, Boundaries &sbc, double energy);

	// ***** construct viscosity ***** //
	static void CalcVis(Flow &vis, const Vctr &vel, double Re, int bftype);

	// ***** construct body forces ***** //
	static void CalcFb(Vctr &fb, const std::vector<double> &mpg);
	static void CalcFb(Vctr &fb, const std::vector<double> &mpg, const Vctr &f);
	static void CalcFb(Vctr &fb, const std::vector<double> &mpg, const char *filename);

	// ***** apply boundary conditions ***** //
	static void SetBoundaries(Vctr &vel, const Boundaries &bc, const Boundaries &sbc);

	// ***** adjusting mean pressure gradient ***** //
	static void CalcMpg(std::vector<double> &mpg, Vctr &vel, double dt);
	static void CalcMpg(std::vector<double> &mpg, Vctr &vel, double dt, const std::vector<double> &mpgref);

	// ***** utilitis ***** //
	double InnerScale() const;

	// ***** data manipulation ***** //
	static void SwapBulk(Vctr &vel);
	static void RollBack(Flow &fld, const Flow &fldh, double dt);
	static void RemoveSpanMean(Vctr &vel);

};




		




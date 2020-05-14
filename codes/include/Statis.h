# pragma once

# include "Basic.h"


class Statis
{
public:
	// Statis() {}; // default constructor, for static function call
	Statis(const Mesh &ms);
	~Statis();

	void Check(const Flow &fld, const Flow &vis, double Re, double dt);
	void WriteProfile(const char *path, int tstep=-1) const;
	void WriteLogfile(const char *path, int tstep, double time, const double mpg[3]) const;
	static double GetLogtime(const char *path, int tstep);
	
private:
	const Mesh &ms;

	double div_; int divpos_[3];
	double cfl_; int cflpos_[3];
	double velm_[3]; // bulk mean velocities in 3 directions
	double taub_[3]; // 3 components of total stress acting on the boundary
	double ener_;    // total fluctuation energy in the domain

	double *um_,  *vm_,  *wm_;
	double *r11_, *r22_, *r33_;
	double *r12_, *r23_, *r13_;
	double *rpu_, *rpv_, *rpw_;
	double *pm_,  *rpp_, *num_;
	
	double CheckProf(const Flow &fld, const Flow &vis, double velm[3]);

	static double CheckDiv(const Vctr &vel, int pos[3]);
	static double CheckCfl(const Vctr &vel, double dt, int pos[3]);
	static void CheckTaub(const Vctr &vel, const Flow &vis, double taub[3]);

	static long int GetLogpos(const char *path, int tstep);
};
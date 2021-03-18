#pragma once

#include "Basic.h"
#include "Bcond.h"


namespace IDM
{
void calc(Flow &fldh, const Flow &fld,
	const Flow &vis, const Vctr &fb, const Boundaries &bc, const Boundaries &sbc, double dt);

// subroutines for computation

// step 1: calculate RHS of momentum equation
void urhs1(Scla &ruh, const Flow &fld, const Flow &vis, const Scla &fbx, const Boundaries &bc, const Boundaries &sbc);
void urhs2(Scla &rvh, const Flow &fld, const Flow &vis, const Scla &fby, const Boundaries &bc, const Boundaries &sbc);
void urhs3(Scla &rwh, const Flow &fld, const Flow &vis, const Scla &fbz, const Boundaries &bc, const Boundaries &sbc);

// step 2: calculate intermedia velocity (solve TDMAs)
void getuh1(Vctr &velh, const Vctr &vel, const Flow &vis, const Boundaries &sbc, double dt);
void getuh2(Vctr &velh, const Vctr &vel, const Flow &vis, const Boundaries &sbc, double dt);
void getuh3(Vctr &velh, const Vctr &vel, const Flow &vis, const Boundaries &sbc, double dt);

// step 3: calculate projector (solve Poisson equation with FFT)
void rhsdp(Scla &rdp, const Vctr &velh, const Boundaries &bc, const Boundaries &sbc, double dt);
void getfdp(Scla &dp, const Boundaries &sbc, double refp);

// step 4: update velocity & pressure fields
void update(Flow &fldh, const Flow &fld, double dt);

} // namespace IDM






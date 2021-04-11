#pragma once

#include "Basic.h"
#include "Bcond.h"


namespace OFW
{

void OffWallVelo(Boundaries &bc, Boundaries &sbc, Vctr &veltmp, const Vctr &vel, const Vctr &veldns, double Ret, double rsclx, double rsclu);
void OffWallVelo(Boundaries &bc, Boundaries &sbc, Vctr &veltmp, const Vctr &vel, double time, double Ret);

void OffWallSubGridUniform    (Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu);
void OffWallSubGridShear      (Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu);
void OffWallSubGridDissipation(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu);

double ReynoldsStressDefect(int j, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu);

} // namespace OFW


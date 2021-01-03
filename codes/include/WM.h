#pragma once

#include "Basic.h"

namespace WM
{

void UniformWallShear(Flow &vis, const Vctr &vel, double tau12);
void LogLawWallShear (Flow &vis, const Vctr &vel, double Ret);

void OffWallSubGridUniform    (Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu);
void OffWallSubGridShear      (Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu);
void OffWallSubGridDissipation(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu);

double ReynoldsStressDefect(int j, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu);

void debug_ShowBoundaryShear(const Vctr &vel, const Flow &vis);
} // namespace WM
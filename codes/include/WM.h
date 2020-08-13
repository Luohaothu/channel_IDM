#pragma once

#include "Basic.h"

namespace WM
{

void UniformWallShear(Flow &vis, const Vctr &vel, double tau12);
void OffWallSGS(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re, double Ret, double rsclx, double rsclu);

void debug_ShowBoundaryShear(const Vctr &vel, const Flow &vis);
} // namespace WM
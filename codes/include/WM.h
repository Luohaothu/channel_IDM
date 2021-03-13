#pragma once

#include "Basic.h"

namespace WM
{

void LogLawWallShear(Flow &vis, const Vctr &vel, double Ret);
void UniformVisShear(Flow &vis, const Vctr &vel, double tau12, double tau23);
void UniformReyShear(Flow &vis, const Vctr &vel, double r12ref, double r23ref, double Re);
void SubGridShear   (Flow &vis, const Vctr &vel, const Vctr &shearsgs, double Re);

void debug_ShowBoundaryShear(const Vctr &vel, const Flow &vis);
} // namespace WM
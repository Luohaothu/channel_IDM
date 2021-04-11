#pragma once

#include "Basic.h"


namespace PIO
{

void PredictBoundary(Vctr &veltmp, const Vctr &vel, const Vctr &velmfu, double Ret, double rsclx, double rsclu);
void PredictBoundary(Vctr &veltmp, const Vctr &vel, double time, double Ret);
void PredictBoundarySGS(Vctr &shearsgs, const Vctr &velout, double Ret);

void ReadSynthe(Vctr &vel, double time, double Ret);

// Vctr Predict(double y, const Vctr &velout, const Vctr &veluni, double Ret, double rsclx, double rsclu, int cnt=-1);

} // namespace PIO


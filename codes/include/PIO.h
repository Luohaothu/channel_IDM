#pragma once

#include "Basic.h"


namespace PIO
{

Vctr BoundaryPredict(const Vctr &vel, const Vctr &velmfu, double Ret, double rsclx, double rsclu);
Vctr Predict(double y, const Vctr &velout, const Vctr &veluni, double Ret, double rsclx, double rsclu, int cnt=-1);

void PredictBoundarySGS(Vctr &shearsgs, const Vctr &velout, double Ret);

} // namespace PIO


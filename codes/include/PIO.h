#pragma once

#include "Basic.h"


namespace PIO
{

Vctr BoundaryPredict(const Vctr &vel, const Vctr &velmfu, double Ret, double rsclx, double rsclu);

} // namespace PIO


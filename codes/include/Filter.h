#pragma once

#include "Basic.h"


namespace Filter
{

double TestFilter(int i, int j, int k, const Scla &q);
double FilterNode(
	double x_, double y_, double z_,
	double dx, double dy, double dz, const Scla &q, int stgtyp);

inline double FilterNodeA(double x, double y, double z, double dx, double dy, double dz, const Scla &q) { return FilterNode(x,y,z,dx,dy,dz,q,0); }
inline double FilterNodeU(double x, double y, double z, double dx, double dy, double dz, const Scla &q) { return FilterNode(x,y,z,dx,dy,dz,q,1); }
inline double FilterNodeV(double x, double y, double z, double dx, double dy, double dz, const Scla &q) { return FilterNode(x,y,z,dx,dy,dz,q,2); }
inline double FilterNodeW(double x, double y, double z, double dx, double dy, double dz, const Scla &q) { return FilterNode(x,y,z,dx,dy,dz,q,3); }

} // namespace Filter
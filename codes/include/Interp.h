#pragma once

#include "Basic.h"


namespace Interp
{
double Shift(double x, double x0, double x1);
double BiSearch(double x, const double *xs, int lo, int hi);

double InterpCell(const Scla &q,
	double x, double xm, double xp, int im, int ip,
	double y, double ym, double yp, int jm, int jp,
	double z, double zm, double zp, int km, int kp );

double InterpNode(double x_, double y_, double z_, const Scla &q, int stgtyp);
void   InterpBulk(Scla &dst, const Scla &src, int stgtyp);

inline double InterpNodeA(double x, double y, double z, const Scla &q) { return InterpNode(x,y,z,q,0); }
inline double InterpNodeU(double x, double y, double z, const Scla &q) { return InterpNode(x,y,z,q,1); }
inline double InterpNodeV(double x, double y, double z, const Scla &q) { return InterpNode(x,y,z,q,2); }
inline double InterpNodeW(double x, double y, double z, const Scla &q) { return InterpNode(x,y,z,q,3); }

inline void InterpBulkA(Scla &dst, const Scla &src) { InterpBulk(dst, src, 0); }
inline void InterpBulkU(Scla &dst, const Scla &src) { InterpBulk(dst, src, 1); }
inline void InterpBulkV(Scla &dst, const Scla &src) { InterpBulk(dst, src, 2); }
inline void InterpBulkW(Scla &dst, const Scla &src) { InterpBulk(dst, src, 3); }

} // namespace Interp



// class Interp
// {
// public:
// 	Interp(const Mesh &ms0, const Mesh &ms1);
// 	~Interp();

// 	void InitIndices();

// 	const Mesh &ms0;
// 	const Mesh &ms1;

// 	int *isu, *isa;
// 	int *jsv, *jsa;
// 	int *ksw, *ksa;
// };



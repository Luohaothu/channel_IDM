#include "Basic.h"

using namespace std;


Vctr::Vctr(const Mesh &ms):
ms(ms),
v1_(ms),
v2_(ms),
v3_(ms)
{}


double Vctr::Module(int i, int j, int k) const
/* module of the vector at cell-center (i,j,k) */
{
	int id =        ms.idx(i,j,k);
	int ip, jp, kp; ms.ipx(i,j,k,ip,jp,kp);

	const Scla &u = v1_;
	const Scla &v = v2_;
	const Scla &w = v3_;

	return	.5 * sqrt(
		    pow(u[id] + u[ip], 2.)
		+   pow(v[id] + v[jp], 2.)
		+   pow(w[id] + w[kp], 2.) );
}

double Vctr::Divergence(int i, int j, int k) const
/* compute divergence of a vector field at the center of cell (i,j,k) */
{
	int id =        ms.idx(i,j,k);
	int ip, jp, kp; ms.ipx(i,j,k,ip,jp,kp);

	const Scla &u = v1_;
	const Scla &v = v2_;
	const Scla &w = v3_;

	return (u[ip] - u[id]) / ms.dx(i)
		 + (v[jp] - v[id]) / ms.dy(j)
		 + (w[kp] - w[id]) / ms.dz(k);
}

double Vctr::Convection(int i, int j, int k) const
/* compute convection coefficient of a vector field at the center of cell (i,j,k) */
{
	int id =        ms.idx(i,j,k);
	int ip, jp, kp; ms.ipx(i,j,k,ip,jp,kp);

	const Scla &u = v1_;
	const Scla &v = v2_;
	const Scla &w = v3_;

	return	.5 * fabs(u[ip] + u[id]) / ms.dx(i)
		+	.5 * fabs(v[jp] + v[id]) / ms.dy(j)
		+	.5 * fabs(w[kp] + w[id]) / ms.dz(k);
}

const double* Vctr::ShearStrain(int i, int j, int k) const
// compute the 3 shear components of strain rate tensor on cell edges up to virtual boundary
{
	int id =              ms.idx(i,j,k);
	int im, jm, km;       ms.imx(i,j,k,im,jm,km);
	double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);

	const Scla &u = v1_;
	const Scla &v = v2_;
	const Scla &w = v3_;

	static double sr[3]; // will be overwritten even called from different objects of this class

	sr[0] = .5 * ((u[id]-u[jm]) / hyc + (v[id]-v[im]) / hxc); // S_12
	sr[1] = .5 * ((v[id]-v[km]) / hzc + (w[id]-w[jm]) / hyc); // S_23
	sr[2] = .5 * ((u[id]-u[km]) / hzc + (w[id]-w[im]) / hxc); // S_13

	return sr;
}

const double* Vctr::Strainrate(int i, int j, int k) const
/* compute the strain rate tensor of a vector field at the center of cell (i,j,k) */
// CAUTION: avoid successive calling to this function, because the static return variable will be overwritten every time
{
	int id =              ms.idx(i,j,k);
	int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
	int im, jm, km;       ms.imx(i,j,k,im,jm,km);
	int ipjp, jpkp, ipkp; ms.ppx(i,j,k,ipjp,jpkp,ipkp);
	int ipjm, jpkm, ipkm; ms.pmx(i,j,k,ipjm,jpkm,ipkm);
	int imjp, jmkp, imkp; ms.mpx(i,j,k,imjp,jmkp,imkp);

	double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
	double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
	double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
	double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
	double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);

	double u1, u2, v1, v2, w1, w2;

	const Scla &u = v1_;
	const Scla &v = v2_;
	const Scla &w = v3_;

	static double sr[6]; // will be overwritten even called from different objects of this class

	// S_ii
	sr[0] = (u[ip] - u[id]) / dxc;
	sr[1] = (v[jp] - v[id]) / dyc;
	sr[2] = (w[kp] - w[id]) / dzc;
	// S_12
	v2 = .25/hxp * ((v[id]+v[jp]) * dxp + (v[ip]+v[ipjp]) * dxc);
	v1 = .25/hxc * ((v[id]+v[jp]) * dxm + (v[im]+v[imjp]) * dxc);
	u2 = .25/hyp * ((u[id]+u[ip]) * dyp + (u[jp]+u[ipjp]) * dyc);
	u1 = .25/hyc * ((u[id]+u[ip]) * dym + (u[jm]+u[ipjm]) * dyc);
	sr[3] = .5 * ((v2-v1) / dxc + (u2-u1) / dyc);
	// S_23
	w2 = .25/hyp * ((w[id]+w[kp]) * dyp + (w[jp]+w[jpkp]) * dyc);
	w1 = .25/hyc * ((w[id]+w[kp]) * dym + (w[jm]+w[jmkp]) * dyc);
	v2 = .25/hzp * ((v[id]+v[jp]) * dzp + (v[kp]+v[jpkp]) * dzc);
	v1 = .25/hzc * ((v[id]+v[jp]) * dzm + (v[km]+v[jpkm]) * dzc);
	sr[4] = .5 * ((w2-w1) / dyc + (v2-v1) / dzc);
	// S_13
	w2 = .25/hxp * ((w[id]+w[kp]) * dxp + (w[ip]+w[ipkp]) * dxc);
	w1 = .25/hxc * ((w[id]+w[kp]) * dxm + (w[im]+w[imkp]) * dxc);
	u2 = .25/hzp * ((u[id]+u[ip]) * dzp + (u[kp]+u[ipkp]) * dzc);
	u1 = .25/hzc * ((u[id]+u[ip]) * dzm + (u[km]+u[ipkm]) * dzc);
	sr[5] = .5 * ((w2-w1) / dxc + (u2-u1) / dzc);

	return sr;
}

const double* Vctr::Gradient(int i, int j, int k) const
/* compute the gradient tensor of a vector field at the center of cell (i,j,k) */
// CAUTION: avoid successive calling to this function, because the static return variable will be overwritten every time
{
	int id =              ms.idx(i,j,k);
	int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
	int im, jm, km;       ms.imx(i,j,k,im,jm,km);
	int ipjp, jpkp, ipkp; ms.ppx(i,j,k,ipjp,jpkp,ipkp);
	int ipjm, jpkm, ipkm; ms.pmx(i,j,k,ipjm,jpkm,ipkm);
	int imjp, jmkp, imkp; ms.mpx(i,j,k,imjp,jmkp,imkp);
	int imjm, jmkm, imkm; ms.mmx(i,j,k,imjm,jmkm,imkm);

	double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
	double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
	double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
	double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
	double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);
	double hxm, hym, hzm; ms.hmx(i,j,k,hxm,hym,hzm);
	
	double u1, u2, v1, v2, w1, w2;
	
	const Scla &u = v1_;
	const Scla &v = v2_;
	const Scla &w = v3_;
	
	static double gr[9]; // will be overwritten even called from different objects of this class

	// a_ii
	gr[0] = (u[ip] - u[id]) / dxc;	// a11
	gr[4] = (v[jp] - v[id]) / dyc;	// a22
	gr[8] = (w[kp] - w[id]) / dzc;	// a33

	// a_1i
	v2 = .25/hxp * ((v[id]+v[jp]) * dxp + (v[ip]+v[ipjp]) * dxc);
	v1 = .25/hxc * ((v[id]+v[jp]) * dxm + (v[im]+v[imjp]) * dxc);
	w2 = .25/hxp * ((w[id]+w[kp]) * dxp + (w[ip]+w[ipkp]) * dxc);
	w1 = .25/hxc * ((w[id]+w[kp]) * dxm + (w[im]+w[imkp]) * dxc);
	gr[1] = (v2-v1) / dxc;	// a12
	gr[2] = (w2-w1) / dxc;	// a13

	// a_2i
	u2 = .25/hyp * ((u[id]+u[ip]) * dyp + (u[jp]+u[ipjp]) * dyc);
	u1 = .25/hyc * ((u[id]+u[ip]) * dym + (u[jm]+u[ipjm]) * dyc);
	w2 = .25/hyp * ((w[id]+w[kp]) * dyp + (w[jp]+w[jpkp]) * dyc);
	w1 = .25/hyc * ((w[id]+w[kp]) * dym + (w[jm]+w[jmkp]) * dyc);
	gr[3] = (u2-u1) / dyc;	// a21
	gr[5] = (w2-w1) / dyc;	// a23
	
	// a_3i
	u2 = .25/hzp * ((u[id]+u[ip]) * dzp + (u[kp]+u[ipkp]) * dzc);
	u1 = .25/hzc * ((u[id]+u[ip]) * dzm + (u[km]+u[ipkm]) * dzc);
	v2 = .25/hzp * ((v[id]+v[jp]) * dzp + (v[kp]+v[jpkp]) * dzc);
	v1 = .25/hzc * ((v[id]+v[jp]) * dzm + (v[km]+v[jpkm]) * dzc);
	gr[6] = (u2-u1) / dzc;	// a31
	gr[7] = (v2-v1) / dzc;	// a32

	return gr;
}





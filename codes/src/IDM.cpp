#include "IDM.h"
#include "Matrix.h"

using namespace std;


void IDM::calc(Flow &fld, Flow &fldh,
	const Flow &vis, const Vctr &fb, const Boundaries &bc, const Boundaries &sbc, double dt)
{
	const Vctr &velh = fldh.GetVec();
	const Vctr &vel = fld.SeeVec();
	const Scla &p = fld.SeeScl();
	Scla &dp = fldh.GetScl();

	// uhcalc
	#pragma omp parallel
	{
		urhs1(fldh.GetVec(1), fld, vis, fb[1], bc, sbc);
		urhs2(fldh.GetVec(2), fld, vis, fb[2], bc, sbc);
		urhs3(fldh.GetVec(3), fld, vis, fb[3], bc, sbc);

		#pragma omp single
		fldh.CleanBoundary();
		
		getuh1(fldh.GetVec(), vel, vis, sbc, dt);
		getuh2(fldh.GetVec(), vel, vis, sbc, dt);
		getuh3(fldh.GetVec(), vel, vis, sbc, dt);
	}

	// dpcalc
	#pragma omp parallel
	rhsdp(dp, velh, bc, sbc, dt); // rdp (which shares memory with dp)
	
	if (fld.ms.x(0)) dp.dctxz(); // rdp->frdp
	else             dp.fftxz();

	#pragma omp parallel
	getfdp(dp, sbc, p.meanxz(fld.ms.Ny-1)); // frdp->fdp

	if (fld.ms.x(0)) dp.idctxz(); // fdp->dp
	else             dp.ifftxz();

	// upcalc
	#pragma omp parallel
	update(fld, fldh, dt);
}


/***** intermediate velocity computation *****/

void IDM::urhs1(Scla &ruh,
	const Flow &fld, const Flow &vis, const Scla &fbx, const Boundaries &bc, const Boundaries &sbc)
/* compute right hand side R_1 for intermediate velocity at all no-boundary grid points */
{
	double u1, u2, v1, v2, w1, w2;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double api, aci, ami, apj, acj, amj, apk, ack, amk;
	double l11un, l12vn, l13wn, m11un, m12vn, m13wn, pressg;

	bool   ifb1, ifb2, ifb3, ifb4;
	double ubc1, ubc2, ubc3, ubc4, vbc3, vbc4;
	double sbu1, sbu2, sbu3, sbu4, sbv3, sbv4;

	const Mesh &ms = ruh.ms;
	const Scla &u = fld.SeeVec(1), &nux = vis.SeeVec(1);
	const Scla &v = fld.SeeVec(2), &nuy = vis.SeeVec(2);
	const Scla &w = fld.SeeVec(3), &nuz = vis.SeeVec(3);
	const Scla &p = fld.SeeScl(),  &nuc = vis.SeeScl();

	#pragma omp for
	for (int j=1; j<ms.Ny; j++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hym,hyp,hyc=ms.hy(j,hym,hyp);
	for (int k=1; k<ms.Nz; k++) { double dzm,dzp,dzc=ms.dz(k,dzm,dzp), hzm,hzp,hzc=ms.hz(k,hzm,hzp);
	for (int i=1; i<ms.Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxm,hxp,hxc=ms.hx(i,hxm,hxp);

		int id =              ms.idx(i,j,k);
		int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
		int im, jm, km;       ms.imx(i,j,k,im,jm,km);
		int imjp, jmkp, imkp; ms.mpx(i,j,k,imjp,jmkp,imkp);

		double mbcx = 0, mbcy = 0; // mbcx is the mbc induced by differencing in direction x (by operator M^1)

		if (ifb1 = (ms.x(0) && i==2)      ) { ubc1 = bc.ub1(j,k); sbu1 = sbc.ub1(j,k); }
		if (ifb2 = (ms.x(0) && i==ms.Nx-1)) { ubc2 = bc.ub2(j,k); sbu2 = sbc.ub2(j,k); }
		if (ifb3 = (ms.y(0) && j==1)      ) { ubc3 = bc.ub3(i,k); sbu3 = sbc.ub3(i,k); }
		if (ifb4 = (ms.y(0) && j==ms.Ny-1)) { ubc4 = bc.ub4(i,k); sbu4 = sbc.ub4(i,k); }

		if (ifb3) { vbc3 = .5/hxc * ( bc.vb3(i,k) * dxm +  bc.vb3(ms.ima(i),k) * dxc);
		            sbv3 = .5/hxc * (sbc.vb3(i,k) * dxm + sbc.vb3(ms.ima(i),k) * dxc); }
		if (ifb4) { vbc4 = .5/hxc * ( bc.vb4(i,k) * dxm +  bc.vb4(ms.ima(i),k) * dxc);
		            sbv4 = .5/hxc * (sbc.vb4(i,k) * dxm + sbc.vb4(ms.ima(i),k) * dxc); }

		// interpolate viscosity & velocity at step n to 6 neighbour points around ruh[id]
		vis1 = nuc[im]; vis2 = nuc[id];
		vis3 = nuz[id]; vis4 = nuz[jp];
		vis5 = nuy[id]; vis6 = nuy[kp];

		u1 = .5 * (u[id] + u[im]);
		u2 = .5 * (u[id] + u[ip]);
		v1 = .5/hxc * (v[id]*dxm + v[im]  *dxc);
		v2 = .5/hxc * (v[jp]*dxm + v[imjp]*dxc);
		w1 = .5/hxc * (w[id]*dxm + w[im]  *dxc);
		w2 = .5/hxc * (w[kp]*dxm + w[imkp]*dxc);

		// viscous terms
		api = 1./hxc * vis2/dxc;
		aci =-1./hxc *(vis2/dxc + vis1/dxm);
		ami = 1./hxc * vis1/dxm;
		apj = .5/dyc * vis4/hyp;
		acj =-.5/dyc *(vis4/hyp + vis3/hyc);
		amj = .5/dyc * vis3/hyc;
		apk = .5/dzc * vis6/hzp;
		ack =-.5/dzc *(vis6/hzp + vis5/hzc);
		amk = .5/dzc * vis5/hzc;

		l11un =	api*u[ip] + aci*u[id] + ami*u[im]
		      + apj*u[jp] + acj*u[id] + amj*u[jm]
		      + apk*u[kp] + ack*u[id] + amk*u[km];
		l12vn = .5/hxc/dyc * (vis4 * (v[jp]-v[imjp]) - vis3 * (v[id]-v[im]));
		l13wn = .5/hxc/dzc * (vis6 * (w[kp]-w[imkp]) - vis5 * (w[id]-w[im]));

		// non-linear terms
		// m11un
		api = .5 /hxc * u2                       - api;
		aci = .5 /hxc *(u2 - u1)                 - aci;
		ami =-.5 /hxc * u1                       - ami;
		apj = .25/hyp * v2                       - apj;
		acj = .25/dyc *(v2*dyp/hyp - v1*dym/hyc) - acj;
		amj =-.25/hyc * v1                       - amj;
		apk = .25/hzp * w2                       - apk;
		ack = .25/dzc *(w2*dzp/hzp - w1*dzm/hzc) - ack;
		amk =-.25/hzc * w1                       - amk;

		if (ifb1 || ifb2 || ifb3 || ifb4) {
			// boundary coefficients
			mbcx-= ifb2 * api * ubc2 + ifb1 * ami * ubc1;
			mbcy-= ifb4 * apj * ubc4 + ifb3 * amj * ubc3;

			aci -= ifb2 * api * sbu2 + ifb1 * ami * sbu1;
			acj -= ifb4 * apj * sbu4 + ifb3 * amj * sbu3;
			api -= ifb2 * api;
			apj -= ifb4 * apj;
			ami -= ifb1 * ami;
			amj -= ifb3 * amj;
		}

		m11un =	api*u[ip] + aci*u[id] + ami*u[im]
		      + apj*u[jp] + acj*u[id] + amj*u[jm]
		      + apk*u[kp] + ack*u[id] + amk*u[km];

		// m12vn
		u2 = .5/hyp * (u[id]*dyp + u[jp]*dyc);
		u1 = .5/hyc * (u[id]*dym + u[jm]*dyc);

		if (ifb3 || ifb4) {
			mbcy -= .5/dyc * (
				ifb4 * (u2*vbc4 - vis4/hxc * (bc.vb4(i,k)-bc.vb4(ms.ima(i),k))) -
				ifb3 * (u1*vbc3 - vis3/hxc * (bc.vb3(i,k)-bc.vb3(ms.ima(i),k))) );

			u1 += ifb4 * u2 * sbv4; vis3 += ifb4 * vis4 * sbv4;
			u2 += ifb3 * u1 * sbv3; vis4 += ifb3 * vis3 * sbv3;
			u2 -= ifb4 * u2;        vis4 -= ifb4 * vis4;
			u1 -= ifb3 * u1;        vis3 -= ifb3 * vis3;
			// put part of M's boundary midification into L: L\~ - M = L\~ - (N - L) = (L + L\~)/2 - [ N - (L + L\~)/2 ]
			l12vn += .5 * (.5/hxc/dyc * (vis4 * (v[jp]-v[imjp]) - vis3 * (v[id]-v[im])) - l12vn);
		}

		m12vn = .5/dyc * (u2*v2 - u1*v1) - l12vn;

		// m13wn
		u2 = .5/hzp * (u[id]*dzp + u[kp]*dzc);
		u1 = .5/hzc * (u[id]*dzm + u[km]*dzc);

		m13wn = .5/dzc * (u2*w2 - u1*w1) - l13wn;

		// pressure gradient term
		pressg = (p[id] - p[im]) / hxc;

		// R_1 with boundary modification
		ruh[id] = (l11un + l12vn + l13wn)
		        - (m11un + m12vn + m13wn)
		        - pressg + fbx[id] + mbcx + mbcy;
	}}}
}

void IDM::urhs2(Scla &rvh,
	const Flow &fld, const Flow &vis, const Scla &fby, const Boundaries &bc, const Boundaries &sbc)
/* compute right hand side R_2 for intermediate velocity at all non-boundary grid points */
{
	double u1, u2, v1, v2, w1, w2;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double api, aci, ami, apj, acj, amj, apk, ack, amk;
	double l21un, l22vn, l23wn, m21un, m22vn, m23wn, pressg;

	bool   ifb1, ifb2, ifb3, ifb4;
	double ubc1, ubc2, vbc1, vbc2, vbc3, vbc4;
	double sbu1, sbu2, sbv1, sbv2, sbv3, sbv4;

	const Mesh &ms = rvh.ms;
	const Scla &u = fld.SeeVec(1), &nux = vis.SeeVec(1);
	const Scla &v = fld.SeeVec(2), &nuy = vis.SeeVec(2);
	const Scla &w = fld.SeeVec(3), &nuz = vis.SeeVec(3);
	const Scla &p = fld.SeeScl(),  &nuc = vis.SeeScl();

	#pragma omp for
	for (int j=1; j<ms.Ny; j++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hym,hyp,hyc=ms.hy(j,hym,hyp);
	for (int k=1; k<ms.Nz; k++) { double dzm,dzp,dzc=ms.dz(k,dzm,dzp), hzm,hzp,hzc=ms.hz(k,hzm,hzp);
	for (int i=1; i<ms.Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxm,hxp,hxc=ms.hx(i,hxm,hxp);

		int id   =      ms.idx(i,j,k);
		int ipjm =      ms.idx(ms.ipa(i),ms.jma(j),k);
		int jmkp =      ms.idx(i,ms.jma(j),ms.kpa(k));
		int ip, jp, kp; ms.ipx(i,j,k,ip,jp,kp);
		int im, jm, km; ms.imx(i,j,k,im,jm,km);

		double mbcx = 0, mbcy = 0;

		if (ifb1 = (ms.x(0) && i==1)      ) { vbc1 = bc.vb1(j,k); sbv1 = sbc.vb1(j,k); }
		if (ifb2 = (ms.x(0) && i==ms.Nx-1)) { vbc2 = bc.vb2(j,k); sbv2 = sbc.vb2(j,k); }
		if (ifb3 = (ms.y(0) && j==2)      ) { vbc3 = bc.vb3(i,k); sbv3 = sbc.vb3(i,k); }
		if (ifb4 = (ms.y(0) && j==ms.Ny-1)) { vbc4 = bc.vb4(i,k); sbv4 = sbc.vb4(i,k); }

		if (ifb1) { ubc1 = .5/hyc * ( bc.ub1(j,k) * dym +  bc.ub1(ms.jma(j),k) * dyc);
		            sbu1 = .5/hyc * (sbc.ub1(j,k) * dym + sbc.ub1(ms.jma(j),k) * dyc); }
		if (ifb2) { ubc2 = .5/hyc * ( bc.ub2(j,k) * dym +  bc.ub2(ms.jma(j),k) * dyc);
		            sbu2 = .5/hyc * (sbc.ub2(j,k) * dym + sbc.ub2(ms.jma(j),k) * dyc); }

		// interpolate viscosity & velocity at step n to 6 neighbour points around rvh[id]
		vis1 = nuz[id]; vis2 = nuz[ip];
		vis3 = nuc[jm]; vis4 = nuc[id];
		vis5 = nux[id]; vis6 = nux[kp];

		u2 = .5/hyc * (u[ip]*dym + u[ipjm]*dyc);
		u1 = .5/hyc * (u[id]*dym + u[jm]  *dyc);
		v2 = .5 * (v[id] + v[jp]);
		v1 = .5 * (v[id] + v[jm]);
		w2 = .5/hyc * (w[kp]*dym + w[jmkp]*dyc);
		w1 = .5/hyc * (w[id]*dym + w[jm]  *dyc);

		// viscous terms
		api = .5/dxc * vis2/hxp;
		aci =-.5/dxc *(vis2/hxp + vis1/hxc);
		ami = .5/dxc * vis1/hxc;
		apj = 1./hyc * vis4/dyc;
		acj =-1./hyc *(vis4/dyc + vis3/dym);
		amj = 1./hyc * vis3/dym;
		apk = .5/dzc * vis6/hzp;
		ack =-.5/dzc *(vis6/hzp + vis5/hzc);
		amk = .5/dzc * vis5/hzc;

		l22vn =	api*v[ip] + aci*v[id] + ami*v[im]
		      + apj*v[jp] + acj*v[id] + amj*v[jm]
		      + apk*v[kp] + ack*v[id] + amk*v[km];
		l21un = .5/hyc/dxc * (vis2 * (u[ip]-u[ipjm]) - vis1 * (u[id]-u[jm]));
		l23wn = .5/hyc/dzc * (vis6 * (w[kp]-w[jmkp]) - vis5 * (w[id]-w[jm]));

		// non-linear terms
		// m22vn
		api = .25/hxp * u2                       - api;
		aci = .25/dxc *(u2*dxp/hxp - u1*dxm/hxc) - aci;
		ami =-.25/hxc * u1                       - ami;
		apj = .5 /hyc * v2                       - apj;
		acj = .5 /hyc *(v2 - v1)                 - acj;
		amj =-.5 /hyc * v1                       - amj;
		apk = .25/hzp * w2                       - apk;
		ack = .25/dzc *(w2*dzp/hzp - w1*dzm/hzc) - ack;
		amk =-.25/hzc * w1                       - amk;

		if (ifb1 || ifb2 || ifb3 || ifb4) {
			mbcx-= ifb2 * api * vbc2 + ifb1 * ami * vbc1;
			mbcy-= ifb4 * apj * vbc4 + ifb3 * amj * vbc3;

			aci -= ifb2 * api * sbv2 + ifb1 * ami * sbv1;
			acj -= ifb4 * apj * sbv4 + ifb3 * amj * sbv3;
			api -= ifb2 * api;
			apj -= ifb4 * apj;
			ami -= ifb1 * ami;
			amj -= ifb3 * amj;
		}

		m22vn =	api*v[ip] + aci*v[id] + ami*v[im]
		      + apj*v[jp] + acj*v[id] + amj*v[jm]
		      + apk*v[kp] + ack*v[id] + amk*v[km];

		// m21un
		v2 = .5/hxp * (v[id]*dxp + v[ip]*dxc);
		v1 = .5/hxc * (v[id]*dxm + v[im]*dxc);

		if (ifb1 || ifb2) {
			mbcx -= .5/dxc * (
				ifb2 * (v2*ubc2 - vis2/hyc * (bc.ub2(j,k)-bc.ub2(ms.jma(j),k))) -
				ifb1 * (v1*ubc1 - vis1/hyc * (bc.ub1(j,k)-bc.ub1(ms.jma(j),k))) );

			v1 += ifb2 * v2 * sbu2; vis1 += ifb2 * vis2 * sbu2;
			v2 += ifb1 * v1 * sbu1; vis2 += ifb1 * vis1 * sbu1;
			v2 -= ifb2 * v2;        vis2 -= ifb2 * vis2;
			v1 -= ifb1 * v1;        vis1 -= ifb1 * vis1;
			// put part of M's boundary midification into L: L\~ - M = L\~ - (N - L) = (L + L\~)/2 - [ N - (L + L\~)/2 ]
			l21un += .5 * (.5/hyc/dxc * (vis2 * (u[ip]-u[ipjm]) - vis1 * (u[id]-u[jm])) - l21un);
		}

		m21un = .5/dxc * (v2*u2 - v1*u1) - l21un;

		// m23wn
		v2 = .5/hzp * (v[id]*dzp + v[kp]*dzc);
		v1 = .5/hzc * (v[id]*dzm + v[km]*dzc);

		m23wn = .5/dzc * (v2*w2 - v1*w1) - l23wn;

		// pressure gradient term
		pressg = (p[id] - p[jm]) / hyc;

		// R_2 with boundary modification
		rvh[id] = (l21un + l22vn + l23wn)
		        - (m21un + m22vn + m23wn)
		        - pressg + fby[id] + mbcx + mbcy;
	}}}
}

void IDM::urhs3(Scla &rwh,
	const Flow &fld, const Flow &vis, const Scla &fbz, const Boundaries &bc, const Boundaries &sbc)
/* compute right hand side R_3 for intermediate velocity at all non-boundary grid points */
{
	double u1, u2, v1, v2, w1, w2;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double api, aci, ami, apj, acj, amj, apk, ack, amk;
	double l31un, l32vn, l33wn, m31un, m32vn, m33wn, pressg;

	bool   ifb1, ifb2, ifb3, ifb4;
	double ubc1, ubc2, vbc3, vbc4, wbc1, wbc2, wbc3, wbc4;
	double sbu1, sbu2, sbv3, sbv4, sbw1, sbw2, sbw3, sbw4;

	const Mesh &ms = rwh.ms;
	const Scla &u = fld.SeeVec(1), &nux = vis.SeeVec(1);
	const Scla &v = fld.SeeVec(2), &nuy = vis.SeeVec(2);
	const Scla &w = fld.SeeVec(3), &nuz = vis.SeeVec(3);
	const Scla &p = fld.SeeScl(),  &nuc = vis.SeeScl();

	#pragma omp for
	for (int j=1; j<ms.Ny; j++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hym,hyp,hyc=ms.hy(j,hym,hyp);
	for (int k=1; k<ms.Nz; k++) { double dzm,dzp,dzc=ms.dz(k,dzm,dzp), hzm,hzp,hzc=ms.hz(k,hzm,hzp);
	for (int i=1; i<ms.Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxm,hxp,hxc=ms.hx(i,hxm,hxp);

		int id =              ms.idx(i,j,k);
		int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
		int im, jm, km;       ms.imx(i,j,k,im,jm,km);
		int ipjm, jpkm, ipkm; ms.pmx(i,j,k,ipjm,jpkm,ipkm);

		double mbcx = 0, mbcy = 0;

		if (ifb1 = (ms.x(0) && i==1)      ) { wbc1 = bc.wb1(j,k); sbw1 = sbc.wb1(j,k); }
		if (ifb2 = (ms.x(0) && i==ms.Nx-1)) { wbc2 = bc.wb2(j,k); sbw2 = sbc.wb2(j,k); }
		if (ifb3 = (ms.y(0) && j==1)      ) { wbc3 = bc.wb3(i,k); sbw3 = sbc.wb3(i,k); }
		if (ifb4 = (ms.y(0) && j==ms.Ny-1)) { wbc4 = bc.wb4(i,k); sbw4 = sbc.wb4(i,k); }

		if (ifb1) { ubc1 = .5/hzc * ( bc.ub1(j,k) * dzm +  bc.ub1(j,ms.kma(k)) * dzc);
		            sbu1 = .5/hzc * (sbc.ub1(j,k) * dzm + sbc.ub1(j,ms.kma(k)) * dzc); }
		if (ifb2) { ubc2 = .5/hzc * ( bc.ub2(j,k) * dzm +  bc.ub2(j,ms.kma(k)) * dzc);
		            sbu2 = .5/hzc * (sbc.ub2(j,k) * dzm + sbc.ub2(j,ms.kma(k)) * dzc); }
		if (ifb3) { vbc3 = .5/hzc * ( bc.vb3(i,k) * dzm +  bc.vb3(i,ms.kma(k)) * dzc);
		            sbv3 = .5/hzc * (sbc.vb3(i,k) * dzm + sbc.vb3(i,ms.kma(k)) * dzc); }
		if (ifb4) { vbc4 = .5/hzc * ( bc.vb4(i,k) * dzm +  bc.vb4(i,ms.kma(k)) * dzc);
		            sbv4 = .5/hzc * (sbc.vb4(i,k) * dzm + sbc.vb4(i,ms.kma(k)) * dzc); }

		// interpolate viscosity and velocity at step n to the position needed
		vis1 = nuy[id]; vis2 = nuy[ip];
		vis3 = nux[id]; vis4 = nux[jp];
		vis5 = nuc[km]; vis6 = nuc[id];

		u2 = .5/hzc * (u[ip]*dzm + u[ipkm]*dzc);
		u1 = .5/hzc * (u[id]*dzm + u[km]  *dzc);
		v2 = .5/hzc * (v[jp]*dzm + v[jpkm]*dzc);
		v1 = .5/hzc * (v[id]*dzm + v[km]  *dzc);
		w2 = .5 * (w[id] + w[kp]);
		w1 = .5 * (w[id] + w[km]);

		// viscous terms
		api = .5/dxc * vis2/hxp;
		aci =-.5/dxc *(vis2/hxp + vis1/hxc);
		ami = .5/dxc * vis1/hxc;
		apj = .5/dyc * vis4/hyp;
		acj =-.5/dyc *(vis4/hyp + vis3/hyc);
		amj = .5/dyc * vis3/hyc;
		apk = 1./hzc * vis6/dzc;
		ack =-1./hzc *(vis6/dzc + vis5/dzm);
		amk = 1./hzc * vis5/dzm;

		l33wn =	api*w[ip] + aci*w[id] + ami*w[im]
		      + apj*w[jp] + acj*w[id] + amj*w[jm]
		      + apk*w[kp] + ack*w[id] + amk*w[km];
		l31un = .5/hzc/dxc * (vis2 * (u[ip]-u[ipkm]) - vis1 * (u[id]-u[km]));
		l32vn = .5/hzc/dyc * (vis4 * (v[jp]-v[jpkm]) - vis3 * (v[id]-v[km]));

		// non-linear terms
		// m33wn
		api = .25/hxp * u2                       - api;
		aci = .25/dxc *(u2*dxp/hxp - u1*dxm/hxc) - aci;
		ami =-.25/hxc * u1                       - ami;
		apj = .25/hyp * v2                       - apj;
		acj = .25/dyc *(v2*dyp/hyp - v1*dym/hyc) - acj;
		amj =-.25/hyc * v1                       - amj;
		apk = .5 /hzc * w2                       - apk;
		ack = .5 /hzc *(w2 - w1)                 - ack;
		amk =-.5 /hzc * w1                       - amk;

		if (ifb1 || ifb2 || ifb3 || ifb4) {
			mbcx-= ifb2 * api * wbc2 + ifb1 * ami * wbc1;
			mbcy-= ifb4 * apj * wbc4 + ifb3 * amj * wbc3;

			aci -= ifb2 * api * sbw2 + ifb1 * ami * sbw1;
			acj -= ifb4 * apj * sbw4 + ifb3 * amj * sbw3;
			api -= ifb2 * api;
			apj -= ifb4 * apj;
			ami -= ifb1 * ami;
			amj -= ifb3 * amj;
		}

		m33wn =	api*w[ip] + aci*w[id] + ami*w[im]
		      + apj*w[jp] + acj*w[id] + amj*w[jm]
		      + apk*w[kp] + ack*w[id] + amk*w[km];

		// m32vn
		w2 = .5/hyp * (w[id]*dyp + w[jp]*dyc);
		w1 = .5/hyc * (w[id]*dym + w[jm]*dyc);

		if (ifb3 || ifb4) {
			mbcy -= .5/dyc * (
				  ifb4 * (w2*vbc4 - vis4/hzc * (bc.vb4(i,k)-bc.vb4(i,ms.kma(k))))
				- ifb3 * (w1*vbc3 - vis3/hzc * (bc.vb3(i,k)-bc.vb3(i,ms.kma(k)))) );

			w1 += ifb4 * w2 * sbv4; vis3 += ifb4 * vis4 * sbv4;
			w2 += ifb3 * w1 * sbv3; vis4 += ifb3 * vis3 * sbv3;
			w2 -= ifb4 * w2;        vis4 -= ifb4 * vis4;
			w1 -= ifb3 * w1;        vis3 -= ifb3 * vis3;
			// put part of M's boundary midification into L: L\~ - M = L\~ - (N - L) = (L + L\~)/2 - [ N - (L + L\~)/2 ]
			l32vn += .5 * (.5/hzc/dyc * (vis4 * (v[jp]-v[jpkm]) - vis3 * (v[id]-v[km])) - l32vn);
		}

		m32vn = .5/dyc * (w2*v2 - w1*v1) - l32vn;

		// m31un
		w2 = .5/hxp * (w[id]*dxp + w[ip]*dxc);
		w1 = .5/hxc * (w[id]*dxm + w[im]*dxc);

		if (ifb1 || ifb2) {
			mbcx -= .5/dxc * (
				ifb2 * (w2*ubc2 - vis2/hzc * (bc.ub2(j,k)-bc.ub2(j,ms.kma(k)))) -
				ifb1 * (w1*ubc1 - vis1/hzc * (bc.ub1(j,k)-bc.ub1(j,ms.kma(k)))) );

			w1 += ifb2 * w2 * sbu2; vis1 += ifb2 * vis2 * sbu2;
			w2 += ifb1 * w1 * sbu1; vis2 += ifb1 * vis1 * sbu1;
			w2 -= ifb2 * w2;        vis2 -= ifb2 * vis2;
			w1 -= ifb1 * w1;        vis1 -= ifb1 * vis1;
			// put part of M's boundary midification into L: L\~ - M = L\~ - (N - L) = (L + L\~)/2 - [ N - (L + L\~)/2 ]
			l31un = .5 * (.5/hzc/dxc * (vis2 * (u[ip]-u[ipkm]) - vis1 * (u[id]-u[km])) - l31un);
		}

		m31un = .5/dxc * (w2*u2 - w1*u1) - l31un;

		// pressure gradient term
		pressg = (p[id] - p[km]) / hzc;

		// R_3 with boundary modification
		rwh[id] = (l31un + l32vn + l33wn)
		        - (m31un + m32vn + m33wn)
		        - pressg + fbz[id] + mbcx + mbcy;
	}}}
}


void IDM::getuh1(Vctr &velh,
	const Vctr &vel, const Flow &vis, const Boundaries &sbc, double dt)
/* compute deltaU^**, result returned by uh (the RHS ruh should be pre stored in uh ) */
{
	double u1, u2, v1, v2, w1, w2;
	double vis1, vis2, vis3, vis4, vis5, vis6;

	bool   ifb1, ifb2, ifb3, ifb4;
	double sbu1, sbu2, sbu3, sbu4;

	const Mesh &ms = velh.ms;
	
	Scla &uh = velh[1];

	const Scla &u = vel[1], &nuc = vis.SeeScl();
	const Scla &v = vel[2], &nuy = vis.SeeVec(2);
	const Scla &w = vel[3], &nuz = vis.SeeVec(3);

	double *ap = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ac = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *am = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ar = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];

	Matrix matx(ms.Nx-1-(bool)ms.x(0));
	Matrix maty(ms.Ny-1);
	Matrix matz(ms.Nz-1);

	// ( I + dt M_11^2 )
	#pragma omp barrier
	#pragma omp for collapse(2)
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double dxm, dxp, dxc = ms.dx(i,dxm,dxp);
		double hxc = ms.hx(i);

		for (int j=1; j<ms.Ny; j++) {

			int id   = ms.idx(i,j,k);
			int im   = ms.idx(ms.ima(i),j,k);
			int jp   = ms.idx(i,ms.jpa(j),k);
			int imjp = ms.idx(ms.ima(i),ms.jpa(j),k);

			double dym, dyp, dyc = ms.dy(j,dym,dyp);
			double hym, hyp, hyc = ms.hy(j,hym,hyp);

			if (ifb3 = (ms.y(0) && j==1)      ) sbu3 = sbc.ub3(i,k);
			if (ifb4 = (ms.y(0) && j==ms.Ny-1)) sbu4 = sbc.ub4(i,k);

			v1 = .5/hxc * (v[id]*dxm + v[im]  *dxc);
			v2 = .5/hxc * (v[jp]*dxm + v[imjp]*dxc);
			vis3 = nuz[id];
			vis4 = nuz[jp];

			ap[j] = ( .25/hyp * v2                       - .5/dyc * vis4/hyp             ) * dt;
			ac[j] = ( .25/dyc *(v2*dyp/hyp - v1*dym/hyc) + .5/dyc *(vis4/hyp + vis3/hyc) ) * dt + 1;
			am[j] = (-.25/hyc * v1                       - .5/dyc * vis3/hyc             ) * dt;

			if (ifb3 || ifb4) ac[j] -= ifb4 * ap[j] * sbu4 + ifb3 * am[j] * sbu3;

			ar[j] = dt * uh[id];
		}

		if (ms.y(0)) maty. tdma(&am[1], &ac[1], &ap[1], &ar[1]); // apj at j=Ny-1 and amj at j=1 are redundant in tdma, thus no need for explicit removal
		else         maty.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int j=1; j<ms.Ny; j++) uh(i,j,k) = ar[j];
	}}

	// ( I + dt M_11^1 )
	#pragma omp barrier
	#pragma omp for collapse(2)
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			int id = ms.idx(i,j,k);
			int im = ms.idx(ms.ima(i),j,k);
			int ip = ms.idx(ms.ipa(i),j,k);

			double dxm, dxp, dxc = ms.dx(i,dxm,dxp);
			double hxm, hxp, hxc = ms.hx(i,hxm,hxp);

			if (ifb1 = (ms.x(0) && i==2)      ) sbu1 = sbc.ub1(j,k);
			if (ifb2 = (ms.x(0) && i==ms.Nx-1)) sbu2 = sbc.ub2(j,k);

			u1 = .5 * (u[id] + u[im]);
			u2 = .5 * (u[id] + u[ip]);
			vis1 = nuc[im];
			vis2 = nuc[id];

			ap[i] = ( .5/hxc * u2       - 1./hxc * vis2/dxc             ) * dt;
			ac[i] = ( .5/hxc *(u2 - u1) + 1./hxc *(vis2/dxc + vis1/dxm) ) * dt + 1;
			am[i] = (-.5/hxc * u1       - 1./hxc * vis1/dxm             ) * dt;

			if (ifb1 || ifb2) ac[i] -= ifb2 * ap[i] * sbu2 + ifb1 * am[i] * sbu1;

			ar[i] = uh[id];
		}

		if (ms.x(0)) matx. tdma(&am[2], &ac[2], &ap[2], &ar[2]);
		else         matx.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int i=1; i<ms.Nx; i++) uh(i,j,k) = ar[i];
	}}

	// ( I + dt M_11^3 )
	#pragma omp barrier
	#pragma omp for collapse(2)
	for (int j=1; j<ms.Ny; j++) {
	for (int i=1; i<ms.Nx; i++) {

		double dxm, dxp, dxc = ms.dx(i,dxm,dxp);
		double hxc = ms.hx(i);

		for (int k=1; k<ms.Nz; k++) {

			int id   = ms.idx(i,j,k);
			int im   = ms.idx(ms.ima(i),j,k);
			int kp   = ms.idx(i,j,ms.kpa(k));
			int imkp = ms.idx(ms.ima(i),j,ms.kpa(k));

			double dzm, dzp, dzc = ms.dz(k,dzm,dzp);
			double hzm, hzp, hzc = ms.hz(k,hzm,hzp);

			w1 = .5/hxc * (w[id]*dxm + w[im]  *dxc);
			w2 = .5/hxc * (w[kp]*dxm + w[imkp]*dxc);
			vis5 = nuy[id];
			vis6 = nuy[kp];

			ap[k] = ( .25/hzp * w2                       - .5/dzc * vis6/hzp             ) * dt;
			ac[k] = ( .25/dzc *(w2*dzp/hzp - w1*dzm/hzc) + .5/dzc *(vis6/hzp + vis5/hzc) ) * dt + 1;
			am[k] = (-.25/hzc * w1                       - .5/dzc * vis5/hzc             ) * dt;

			ar[k] = uh[id];
		}

		if (ms.z(0)) matz. tdma(&am[1], &ac[1], &ap[1], &ar[1]);
		else         matz.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int k=1; k<ms.Nz; k++) uh(i,j,k) = ar[k];
	}}

	delete[] ap;
	delete[] ac;
	delete[] am;
	delete[] ar;
}


void IDM::getuh2(Vctr &velh,
	const Vctr &vel, const Flow &vis, const Boundaries &sbc, double dt)
/* compute deltaU^**, result returned by vh (the RHS rvh should be pre stored in uh ) */
{
	double u1, u2, v1, v2, w1, w2, l21uh, m21uh;
	double vis1, vis2, vis3, vis4, vis5, vis6;

	bool   ifb1, ifb2, ifb3, ifb4;
	double sbu1, sbu2, sbv1, sbv2, sbv3, sbv4;

	const Mesh &ms = velh.ms;
	const Scla &uh = velh[1];
	
	Scla &vh = velh[2];

	const Scla &u = vel[1], &nux = vis.SeeVec(1);
	const Scla &v = vel[2], &nuc = vis.SeeScl();
	const Scla &w = vel[3], &nuz = vis.SeeVec(3);

	double *ap = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ac = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *am = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ar = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];

	Matrix matx(ms.Nx-1);
	Matrix maty(ms.Ny-1-(bool)ms.y(0));
	Matrix matz(ms.Nz-1);

	// ( I + dt M_22^2 )
	#pragma omp barrier
	#pragma omp for collapse(2)
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double dxm, dxp, dxc = ms.dx(i,dxm,dxp);
		double hxm, hxp, hxc = ms.hx(i,hxm,hxp);

		for (int j=1; j<ms.Ny; j++) {

			int id   =      ms.idx(i,j,k);
			int ipjm =      ms.idx(ms.ipa(i),ms.jma(j),k);
			int ip, jp, kp; ms.ipx(i,j,k,ip,jp,kp);
			int im, jm, km; ms.imx(i,j,k,im,jm,km);

			double dym, dyp, dyc = ms.dy(j,dym,dyp);
			double hyc = ms.hy(j);

			if (ifb1 = (ms.x(0) && i==1)      ) sbu1 = .5/hyc * (sbc.ub1(j,k) * dym + sbc.ub1(ms.jma(j),k) * dyc);
			if (ifb2 = (ms.x(0) && i==ms.Nx-1)) sbu2 = .5/hyc * (sbc.ub2(j,k) * dym + sbc.ub2(ms.jma(j),k) * dyc);
			if (ifb3 = (ms.y(0) && j==2)      ) sbv3 = sbc.vb3(i,k);
			if (ifb4 = (ms.y(0) && j==ms.Ny-1)) sbv4 = sbc.vb4(i,k);
		
			v2 = .5 * (v[id] + v[jp]);
			v1 = .5 * (v[id] + v[jm]);
			vis3 = nuc[jm];
			vis4 = nuc[id];

			ap[j] = ( .5/hyc* v2       - 1./hyc * vis4/dyc             ) * dt;
			ac[j] = ( .5/hyc*(v2 - v1) + 1./hyc *(vis4/dyc + vis3/dym) ) * dt + 1;
			am[j] = (-.5/hyc* v1       - 1./hyc * vis3/dym             ) * dt;

			if (ifb3 || ifb4) ac[j] -= ifb4 * ap[j] * sbv4 + ifb3 * am[j] * sbv3;

			// m21uh
			v2 = .5/hxp * (v[id]*dxp + v[ip]*dxc);
			v1 = .5/hxc * (v[id]*dxm + v[im]*dxc);
			vis1 = nuz[id];
			vis2 = nuz[ip];

			if (ifb1 || ifb2) {
				v1 += ifb2 * v2 * sbu2; vis1 += ifb2 * vis2 * sbu2;
				v2 += ifb1 * v1 * sbu1; vis2 += ifb1 * vis1 * sbu1;
				v1 -= ifb1 * v1;        vis1 -= ifb1 * vis1;
				v2 -= ifb2 * v2;        vis2 -= ifb2 * vis2;
			}

			u2 = .5/hyc * (uh[ip]*dym + uh[ipjm]*dyc);
			u1 = .5/hyc * (uh[id]*dym + uh[jm]  *dyc);

			l21uh = .5/hyc/dxc * (vis2 * (uh[ip]-uh[ipjm]) - vis1 * (uh[id]-uh[jm]));
			m21uh = .5/dxc * (v2*u2 - v1*u1) - l21uh;

			ar[j] = dt * (vh[id] - m21uh);
		}

		if (ms.y(0)) maty. tdma(&am[2], &ac[2], &ap[2], &ar[2]); // apj at j=Ny-1 and amj at j=1 are redundant in tdma, thus no need for explicit removal
		else         maty.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);
		
		for (int j=1; j<ms.Ny; j++) vh(i,j,k) = ar[j];
	}}

	// ( I + dt M_22^1 )
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hyc=ms.hy(j);
	for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			int id   = ms.idx(i,j,k);
			int ip   = ms.idx(ms.ipa(i),j,k);
			int jm   = ms.idx(i,ms.jma(j),k);
			int ipjm = ms.idx(ms.ipa(i),ms.jma(j),k);

			double dxm, dxp, dxc = ms.dx(i,dxm,dxp);
			double hxm, hxp, hxc = ms.hx(i,hxm,hxp);

			if (ifb1 = (ms.x(0) && i==1)      ) sbv1 = sbc.vb1(j,k);
			if (ifb2 = (ms.x(0) && i==ms.Nx-1)) sbv2 = sbc.vb2(j,k);

			u2 = .5/hyc * (u[ip]*dym + u[ipjm]*dyc);
			u1 = .5/hyc * (u[id]*dym + u[jm]  *dyc);
			vis1 = nuz[id];
			vis2 = nuz[ip];

			ap[i] = ( .25/hxp * u2                       - .5/dxc * vis2/hxp             ) * dt;
			ac[i] = ( .25/dxc *(u2*dxp/hxp - u1*dxm/hxc) + .5/dxc *(vis2/hxp + vis1/hxc) ) * dt + 1;
			am[i] = (-.25/hxc * u1                       - .5/dxc * vis1/hxc             ) * dt;

			if (ifb1 || ifb2) ac[i] -= ifb2 * ap[i] * sbv2 + ifb1 * am[i] * sbv1;

			ar[i] = vh[id];
		}

		if (ms.x(0)) matx. tdma(&am[1], &ac[1], &ap[1], &ar[1]);
		else         matx.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int i=1; i<ms.Nx; i++) vh(i,j,k) = ar[i];
	}}

	// ( I + dt M_22^3 )
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hyc=ms.hy(j);
	for (int i=1; i<ms.Nx; i++) {
		for (int k=1; k<ms.Nz; k++) {

			int id   = ms.idx(i,j,k);
			int jm   = ms.idx(i,ms.jma(j),k);
			int kp   = ms.idx(i,j,ms.kpa(k));
			int jmkp = ms.idx(i,ms.jma(j),ms.kpa(k));

			double dzm, dzp, dzc = ms.dz(k,dzm,dzp);
			double hzm, hzp, hzc = ms.hz(k,hzm,hzp);

			w2 = .5/hyc * (w[kp]*dym + w[jmkp]*dyc);
			w1 = .5/hyc * (w[id]*dym + w[jm]  *dyc);
			vis5 = nux[id];
			vis6 = nux[kp];

			ap[k] = ( .25/hzp * w2                       - .5/dzc * vis6/hzp             ) * dt;
			ac[k] = ( .25/dzc *(w2*dzp/hzp - w1*dzm/hzc) + .5/dzc *(vis6/hzp + vis5/hzc) ) * dt + 1;
			am[k] = (-.25/hzc * w1                       - .5/dzc * vis5/hzc             ) * dt;

			ar[k] = vh[id];
		}

		if (ms.z(0)) matz. tdma(&am[1], &ac[1], &ap[1], &ar[1]);
		else         matz.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int k=1; k<ms.Nz; k++) vh(i,j,k) = ar[k];
	}}
	
	delete[] ap;
	delete[] ac;
	delete[] am;
	delete[] ar;
}


void IDM::getuh3(Vctr &velh,
	const Vctr &vel, const Flow &vis, const Boundaries &sbc, double dt)
/* compute deltaW^*, deltaV^*, deltaU^*, and update U^*, V^*, W^* (the RHS rwh should be pre stored in wh ) */
{
	double u1, u2, v1, v2, w1, w2;
	double l31uh, l32vh, l23wh, l12vh, l13wh;
	double m31uh, m32vh, m23wh, m12vh, m13wh;
	double vis1, vis2, vis3, vis4, vis5, vis6;

	bool   ifb1, ifb2, ifb3, ifb4;
	double sbu1, sbu2, sbv3, sbv4, sbw1, sbw2, sbw3, sbw4;

	const Mesh &ms = velh.ms;

	Scla &uh = velh[1];
	Scla &vh = velh[2];
	Scla &wh = velh[3];

	const Scla &u = vel[1], &nux = vis.SeeVec(1);
	const Scla &v = vel[2], &nuy = vis.SeeVec(2);
	const Scla &w = vel[3], &nuz = vis.SeeVec(3);
	const Scla              &nuc = vis.SeeScl();

	double *ap = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ac = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *am = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];
	double *ar = new double[max(max(ms.Nx, ms.Ny), ms.Nz)];

	Matrix matx(ms.Nx-1);
	Matrix maty(ms.Ny-1);
	Matrix matz(ms.Nz-1-(bool)ms.z(0));

	// ( I + dt M_33^2 )
	#pragma omp barrier
	#pragma omp for
	for (int k=1; k<ms.Nz; k++) { double dzm,dzp,dzc=ms.dz(k,dzm,dzp), hzc=ms.hz(k);
	for (int i=1; i<ms.Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxm,hxp,hxc=ms.hx(i,hxm,hxp);
		for (int j=1; j<ms.Ny; j++) {

			int id =              ms.idx(i,j,k);
			int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
			int im, jm, km;       ms.imx(i,j,k,im,jm,km);
			int ipjm, jpkm, ipkm; ms.pmx(i,j,k,ipjm,jpkm,ipkm);
			int imjm, jmkm, imkm; ms.mmx(i,j,k,imjm,jmkm,imkm);

			double dym, dyp, dyc = ms.dy(j,dym,dyp);
			double hym, hyp, hyc = ms.hy(j,hym,hyp);

			if (ifb1 = (ms.x(0) && i==1)      ) { sbu1 = .5/hzc * (sbc.ub1(j,k) * dzm + sbc.ub1(j,ms.kma(k)) * dzc); }
			if (ifb2 = (ms.x(0) && i==ms.Nx-1)) { sbu2 = .5/hzc * (sbc.ub2(j,k) * dzm + sbc.ub2(j,ms.kma(k)) * dzc); }
			if (ifb3 = (ms.y(0) && j==1)      ) { sbv3 = sbc.vb3(i,k); sbw3 = sbc.wb3(i,k); }
			if (ifb4 = (ms.y(0) && j==ms.Ny-1)) { sbv4 = sbc.vb4(i,k); sbw4 = sbc.wb4(i,k); }

			v2 = .5/hzc * (v[jp]*dzm + v[jpkm]*dzc);
			v1 = .5/hzc * (v[id]*dzm + v[km]  *dzc);
			vis3 = nux[id];
			vis4 = nux[jp];

			ap[j] = ( .25/hyp * v2                       - .5/dyc * vis4/hyp             ) * dt;
			ac[j] = ( .25/dyc *(v2*dyp/hyp - v1*dym/hyc) + .5/dyc *(vis4/hyp + vis3/hyc) ) * dt + 1;
			am[j] = (-.25/hyc * v1                       - .5/dyc * vis3/hyc             ) * dt;

			if (ifb3 || ifb4) ac[j] -= ifb4 * ap[j] * sbw4 + ifb3 * am[j] * sbw3;

			// m31uh
			w2 = .5/hxp * (w[id]*dxp + w[ip]*dxc);
			w1 = .5/hxc * (w[id]*dxm + w[im]*dxc);
			vis1 = nuy[id];
			vis2 = nuy[ip];

			if (ifb1 || ifb2) {
				w1 += ifb2 * w2 * sbu2; vis1 += ifb2 * vis2 * sbu2;
				w2 += ifb1 * w1 * sbu1; vis2 += ifb1 * vis1 * sbu1;
				w2 -= ifb2 * w2;        vis2 -= ifb2 * vis2;
				w1 -= ifb1 * w1;        vis1 -= ifb1 * vis1;
			}

			u2 = .5/hzc * (uh[ip]*dzm + uh[ipkm]*dzc);
			u1 = .5/hzc * (uh[id]*dzm + uh[km]  *dzc);

			l31uh = .5/hzc/dxc * (vis2 * (uh[ip]-uh[ipkm]) - vis1 * (uh[id]-uh[km]));
			m31uh = .5/dxc * (w2*u2 - w1*u1) - l31uh;

			// m32vh
			w2 = .5/hyp * (w[id]*dyp + w[jp]*dyc);
			w1 = .5/hyc * (w[id]*dym + w[jm]*dyc);

			if (ifb3 || ifb4) {
				w1 += ifb4 * w2 * sbv4; vis3 += ifb4 * vis4 * sbv4;
				w2 += ifb3 * w1 * sbv3; vis4 += ifb3 * vis3 * sbv3;
				w2 -= ifb4 * w2;        vis4 -= ifb4 * vis4;
				w1 -= ifb3 * w1;        vis3 -= ifb3 * vis3;
			}

			v2 = .5/hzc * (vh[jp]*dzm + vh[jpkm]*dzc);
			v1 = .5/hzc * (vh[id]*dzm + vh[km]  *dzc);

			l32vh = .5/hzc/dyc * (vis4 * (vh[jp]-vh[jpkm]) - vis3 * (vh[id]-vh[km]));
			m32vh = .5/dyc * (w2*v2 - w1*v1) - l32vh;

			ar[j] = dt * (wh[id] - m31uh - m32vh);
		}

		if (ms.y(0)) maty. tdma(&am[1], &ac[1], &ap[1], &ar[1]);
		else         maty.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int j=1; j<ms.Ny; j++) wh(i,j,k) = ar[j];
	}}

	// ( I + dt M_33^1 )
	#pragma omp barrier
	#pragma omp for collapse(2)
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {

		double dzm, dzp, dzc = ms.dz(k,dzm,dzp);
		double hzc = ms.hz(k);

		for (int i=1; i<ms.Nx; i++) {

			int id   = ms.idx(i,j,k);
			int ip   = ms.idx(ms.ipa(i),j,k);
			int km   = ms.idx(i,j,ms.kma(k));
			int ipkm = ms.idx(ms.ipa(i),j,ms.kma(k));

			double dxm, dxp, dxc = ms.dx(i,dxm,dxp);
			double hxm, hxp, hxc = ms.hx(i,hxm,hxp);

			if (ifb1 = (ms.x(0) && i==1)      ) sbw1 = sbc.wb1(j,k);
			if (ifb2 = (ms.x(0) && i==ms.Nx-1)) sbw2 = sbc.wb2(j,k);

			u2 = .5/hzc * (u[ip]*dzm + u[ipkm]*dzc);
			u1 = .5/hzc * (u[id]*dzm + u[km]  *dzc);
			vis1 = nuy[id];
			vis2 = nuy[ip];

			ap[i] = ( .25/hxp * u2                       - .5/dxc * vis2/hxp             ) * dt;
			ac[i] = ( .25/dxc *(u2*dxp/hxp - u1*dxm/hxc) + .5/dxc *(vis2/hxp + vis1/hxc) ) * dt + 1;
			am[i] = (-.25/hxc * u1                       - .5/dxc * vis1/hxc             ) * dt;

			if (ifb1 || ifb2) ac[i] -= ifb2 * ap[i] * sbw2 + ifb1 * am[i] * sbw1;

			ar[i] = wh[id];
		}

		if (ms.x(0)) matx. tdma(&am[1], &ac[1], &ap[1], &ar[1]);
		else         matx.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int i=1; i<ms.Nx; i++) wh(i,j,k) = ar[i];
	}}

	// ( I + dt M_33^3 )
	#pragma omp barrier
	#pragma omp for collapse(2)
	for (int j=1; j<ms.Ny; j++) {
	for (int i=1; i<ms.Nx; i++) {
		for (int k=1; k<ms.Nz; k++) {

			int id = ms.idx(i,j,k);
			int kp = ms.idx(i,j,ms.kpa(k));
			int km = ms.idx(i,j,ms.kma(k));

			double dzm, dzp, dzc = ms.dz(k,dzm,dzp);
			double hzc = ms.hz(k);

			w2 = .5 * (w[id] + w[kp]);
			w1 = .5 * (w[id] + w[km]);
			vis5 = nuc[km];
			vis6 = nuc[id];

			ap[k] = ( .5/hzc * w2       - 1./hzc * vis6/dzc             ) * dt;
			ac[k] = ( .5/hzc *(w2 - w1) + 1./hzc *(vis6/dzc + vis5/dzm) ) * dt + 1;
			am[k] = (-.5/hzc * w1       - 1./hzc * vis5/dzm             ) * dt;

			ar[k] = wh[id];
		}

		if (ms.z(0)) matz. tdma(&am[2], &ac[2], &ap[2], &ar[2]);
		else         matz.ctdma(&am[1], &ac[1], &ap[1], &ar[1]);

		for (int k=1; k<ms.Nz; k++) wh(i,j,k) = ar[k];
	}}

	delete[] ap;
	delete[] ac;
	delete[] am;
	delete[] ar;

	// update dvh
	#pragma omp barrier
	#pragma omp for
	for (int j=2; j<ms.Ny; j++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hyc=ms.hy(j);
	for (int k=1; k<ms.Nz; k++) { double dzm,dzp,dzc=ms.dz(k,dzm,dzp), hzm,hzp,hzc=ms.hz(k,hzm,hzp);
	for (int i=1; i<ms.Nx; i++) {

		int id   = ms.idx(i,j,k);
		int jm   = ms.idx(i,ms.jma(j),k);
		int km   = ms.idx(i,j,ms.kma(k));
		int kp   = ms.idx(i,j,ms.kpa(k));
		int jmkp = ms.idx(i,ms.jma(j),ms.kpa(k));

		// m23wh
		v2 = .5/hzp * (v[id]*dzp + v[kp]*dzc);
		v1 = .5/hzc * (v[id]*dzm + v[km]*dzc);
		vis5 = nux[id];
		vis6 = nux[kp];

		w2 = .5/hyc * (wh[kp]*dym + wh[jmkp]*dyc);
		w1 = .5/hyc * (wh[id]*dym + wh[jm]  *dyc);

		l23wh = .5/hyc/dzc * (vis6 * (wh[kp]-wh[jmkp]) - vis5 * (wh[id]-wh[jm]));
		m23wh = .5/dzc * (v2*w2 - v1*w1) - l23wh;

		vh[id] -= dt * m23wh;
	}}}
	
	// update duh
	#pragma omp barrier
	#pragma omp for
	for (int j=1; j<ms.Ny; j++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hym,hyp,hyc=ms.hy(j,hym,hyp);
	for (int k=1; k<ms.Nz; k++) { double dzm,dzp,dzc=ms.dz(k,dzm,dzp), hzm,hzp,hzc=ms.hz(k,hzm,hzp);
	for (int i=1; i<ms.Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxc=ms.hx(i);

		int id =              ms.idx(i,j,k);
		int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
		int im, jm, km;       ms.imx(i,j,k,im,jm,km);
		int imjp, jmkp, imkp; ms.mpx(i,j,k,imjp,jmkp,imkp);
		
		if (ifb3 = (ms.y(0) && j==1)      ) sbv3 = .5/hxc * (sbc.vb3(i,k) * dxm + sbc.vb3(ms.ima(i),k) * dxc);
		if (ifb4 = (ms.y(0) && j==ms.Ny-1)) sbv4 = .5/hxc * (sbc.vb4(i,k) * dxm + sbc.vb4(ms.ima(i),k) * dxc);

		// m12vh
		u2 = .5/hyp * (u[id]*dyp + u[jp]*dyc);
		u1 = .5/hyc * (u[id]*dym + u[jm]*dyc);
		vis3 = nuz[id];
		vis4 = nuz[jp];

		if (ifb3 || ifb4) {
			u1 += ifb4 * u2 * sbv4; vis3 += ifb4 * vis4 * sbv4;
			u2 += ifb3 * u1 * sbv3; vis4 += ifb3 * vis3 * sbv3;
			u1 -= ifb3 * u1;        vis3 -= ifb3 * vis3;
			u2 -= ifb4 * u2;        vis4 -= ifb4 * vis4;
		}

		v1 = .5/hxc * (vh[id]*dxm + vh[im]  *dxc);
		v2 = .5/hxc * (vh[jp]*dxm + vh[imjp]*dxc);

		l12vh = .5/hxc/dyc * (vis4 * (vh[jp]-vh[imjp]) - vis3 * (vh[id]-vh[im]));
		m12vh = .5/dyc * (u2*v2 - u1*v1) - l12vh;

		// m13wh
		u1 = .5/hzc * (u[id]*dzm + u[km]*dzc);
		u2 = .5/hzp * (u[id]*dzp + u[kp]*dzc);
		vis5 = nuy[id];
		vis6 = nuy[kp];

		w1 = .5/hxc * (wh[id]*dxm + wh[im]  *dxc);
		w2 = .5/hxc * (wh[kp]*dxm + wh[imkp]*dxc);

		l13wh = .5/hxc/dzc * (vis6 * (wh[kp]-wh[imkp]) - vis5 * (wh[id]-wh[im]));
		m13wh = .5/dzc * (u2*w2 - u1*w1) - l13wh;

		uh[id] -= dt * (m12vh + m13wh);
	}}}

	// update intermediate velocity field
	#pragma omp barrier
	#pragma omp for collapse(2)
	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		int id = ms.idx(i,j,k);
		uh[id] += u[id];
		vh[id] += v[id];
		wh[id] += w[id];
	}}}
}



/***** projector computation *****/

void IDM::rhsdp(Scla &rdp,
	const Vctr &velh, const Boundaries &bc, const Boundaries &sbc, double dt)
/* compute RHS of Poisson equation and store in dp */
{
	const Mesh &ms = rdp.ms;
	const Scla &uh = velh[1];
	const Scla &vh = velh[2];
	const Scla &wh = velh[3];

	bool   ifb1, ifb2, ifb3, ifb4;
	double ubc1, ubc2, vbc3, vbc4;
	double sbu1, sbu2, sbv3, sbv4;

	#pragma omp for
	for (int j=1; j<ms.Ny; j++) { double dyc=ms.dy(j);
	for (int k=1; k<ms.Nz; k++) { double dzc=ms.dz(k);
	for (int i=1; i<ms.Nx; i++) { double dxc=ms.dx(i);

		int id =        ms.idx(i,j,k);
		int ip, jp, kp; ms.ipx(i,j,k,ip,jp,kp);

		if (ifb1 = (ms.x(0) && i==1)      ) { ubc1 = bc.ub1(j,k); sbu1 = sbc.ub1(j,k); }
		if (ifb2 = (ms.x(0) && i==ms.Nx-1)) { ubc2 = bc.ub2(j,k); sbu2 = sbc.ub2(j,k); }
		if (ifb3 = (ms.y(0) && j==1)      ) { vbc3 = bc.vb3(i,k); sbv3 = sbc.vb3(i,k); }
		if (ifb4 = (ms.y(0) && j==ms.Ny-1)) { vbc4 = bc.vb4(i,k); sbv4 = sbc.vb4(i,k); }

		// ( Du^h - cbc ) / dt
		rdp[id] = 1./dt * (
			1./dxc * ( uh[ip] * (1 + ifb1*sbu1 - ifb2) + ifb2*ubc2 
			         - uh[id] * (1 + ifb2*sbu2 - ifb1) - ifb1*ubc1 )
		+	1./dyc * ( vh[jp] * (1 + ifb3*sbv3 - ifb4) + ifb4*vbc4
			         - vh[id] * (1 + ifb4*sbv4 - ifb3) - ifb3*vbc3 )
		+	1./dzc * ( wh[kp]-wh[id] ) );
	}}}
}

void IDM::getfdp(Scla &fdp, const Boundaries &sbc, double refp)
/* compute FDP (in Fourier space), result returned by fdp (the RHS frdp should be pre stored in fdp ) */
{
	const Mesh &ms = fdp.ms;

	double *ap = new double[ms.Ny];
	double *ac = new double[ms.Ny];
	double *am = new double[ms.Ny];
	double *ar = new double[ms.Ny];
	double *ai = new double[ms.Ny];

	Matrix maty(ms.Ny-1);

	int i0, in, k0, kn;
	// choose wavenumber range to solve by periodicity in x and/or z directions
	if ( ms.x(0) && !ms.z(0)) { i0 = 1; in = ms.Nx; k0 = 0; kn = ms.Nzc; }
	if (!ms.x(0) && !ms.z(0)) { i0 = 0; in = ms.Nxc;k0 = 0; kn = ms.Nz-1; }

	#pragma omp for collapse(2)
	for (int k=k0; k<kn; k++) {
	for (int i=i0; i<in; i++) {

		double sbv3 = sbc.vb3(i,k);
		double sbv4 = sbc.vb4(i,k);

		for (int j=1; j<ms.Ny; j++) {

			double dyc = ms.dy(j);
			double hym, hyp, hyc = ms.hy(j,hym,hyp);

			bool ifb3 = (ms.y(0) && j==1);
			bool ifb4 = (ms.y(0) && j==ms.Ny-1);

			ap[j] = 1./dyc *   (1 + ifb3*sbv3 - ifb4) / hyp;
			ac[j] =-1./dyc * ( (1 + ifb3*sbv3 - ifb4) / hyp
				             + (1 + ifb4*sbv4 - ifb3) / hyc ) - ms.kx2(i) - ms.kz2(k);
			am[j] = 1./dyc *   (1 + ifb4*sbv4 - ifb3) / hyc;

			ar[j] = fdp[ms.idfr(i,j,k)];
			ai[j] = fdp[ms.idfi(i,j,k)];
		}

		// set reference pressure P(kx=0,kz=0,j=1) to 0
		// fdp(kx=0,kz=0,j=1) = -Nxz * fp(kx=0,kz=0,j=1) ==> <dp(j=1)> = - <p(j=1)>
		if (ms.kx2(i) + ms.kz2(k) < INFTSM) {
			ap[ms.Ny-1] = 0;
			ac[ms.Ny-1] = 1;
			am[ms.Ny-1] = 0;
			ar[ms.Ny-1] = - (ms.Nx-1)*(ms.Nz-1) * refp;
			ai[ms.Ny-1] = 0;
		}

		maty.tdma(&am[1], &ac[1], &ap[1], &ar[1]);
		maty.tdma(&am[1], &ac[1], &ap[1], &ai[1]);

		for (int j=1; j<ms.Ny; j++) {
			fdp[ms.idfr(i,j,k)] = ar[j];
			fdp[ms.idfi(i,j,k)] = ai[j];
		}
	}}
	/* note:
		In Huang's program, imag part of FDP at wavenumbers
		(kx=0,kz=0), (kx=0,kz=Nz/2), (kx=Nx/2,kz=0), (kx=Nx/2,kz=Nz/2)
		are explicitly set to 0. This is not needed here, since
		zero imag FRDP lead directly to zero imag FDP	*/

	delete[] ap;
	delete[] ac;
	delete[] am;
	delete[] ar;
	delete[] ai;
}



void IDM::update(Flow &fld, Flow &fldh, double dt)
/* project from U^* to U^n+1 using DP */
{
	const Mesh &ms = fld.ms;

	Scla &u = fld.GetVec(1), &uh = fldh.GetVec(1);
	Scla &v = fld.GetVec(2), &vh = fldh.GetVec(2);
	Scla &w = fld.GetVec(3), &wh = fldh.GetVec(3);
	Scla &p = fld.GetScl(),  &dp = fldh.GetScl();

	#pragma omp for
	for (int j=1; j<ms.Ny; j++) { double hyc=ms.hy(j);
	for (int k=1; k<ms.Nz; k++) { double hzc=ms.hz(k);
	for (int i=1; i<ms.Nx; i++) { double hxc=ms.hx(i);

		int id =        ms.idx(i,j,k);
		int im, jm, km; ms.imx(i,j,k,im,jm,km);

		// store time derivative in intermediate variables
		uh[id] = (uh[id]-u[id]) / dt - (dp[id]-dp[im]) / hxc;
		vh[id] = (vh[id]-v[id]) / dt - (dp[id]-dp[jm]) / hyc;
		wh[id] = (wh[id]-w[id]) / dt - (dp[id]-dp[km]) / hzc;
		
		// update fields
		if (i > (bool)ms.x(0)) u[id] += uh[id] * dt;
		if (j > (bool)ms.y(0)) v[id] += vh[id] * dt;
		if (k > (bool)ms.z(0)) w[id] += wh[id] * dt;
	}}}

	#pragma omp barrier
	#pragma omp for collapse(2)
	for (int j=1; j<ms.Ny; j++)
	for (int k=1; k<ms.Nz; k++)
	for (int i=1; i<ms.Nx; i++)
		p(i,j,k) += dt * (dp(i,j,k) /= dt);
}






//***** codes for time test *****//

// #include <sys/time.h>
// int cnt = 0;
// struct timeval time0, time1;
// double t_urhs=0, t_getuh=0, t_getdp=0, t_update=0;

// void IDM::calc(Flow &fld, Flow &fldh,
// 	const Flow &vis, const Vctr &fb, const Boundaries &bc, const Boundaries &sbc, double dt)
// {
// 	const Vctr &velh = fldh.GetVec();
// 	const Vctr &vel = fld.SeeVec();
// 	const Scla &p = fld.SeeScl();
// 	Scla &dp = fldh.GetScl();

// 	cout << ++ cnt << endl;

// 	gettimeofday(&time0, NULL);
// 		// uhcalc
// 		#pragma omp parallel
// 		{
// 			urhs1(fldh.GetVec(1), fld, vis, fb[1], bc, sbc);
// 			urhs2(fldh.GetVec(2), fld, vis, fb[2], bc, sbc);
// 			urhs3(fldh.GetVec(3), fld, vis, fb[3], bc, sbc);
// 		}
// 	gettimeofday(&time1, NULL); t_urhs += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);

// 	gettimeofday(&time0, NULL);
// 		#pragma omp parallel
// 		{
// 			getuh1(fldh.GetVec(), vel, vis, sbc, dt);
// 			getuh2(fldh.GetVec(), vel, vis, sbc, dt);
// 			getuh3(fldh.GetVec(), vel, vis, sbc, dt);
// 		}
// 	gettimeofday(&time1, NULL); t_getuh += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);

// 	gettimeofday(&time0, NULL);
// 		// dpcalc
// 		#pragma omp parallel
// 		rhsdp(dp, velh, bc, sbc, dt); // rdp (which shares memory with dp)
// 		dp.fftxz();                   // rdp->frdp
// 		#pragma omp parallel
// 		getfdp(dp, sbc, p.meanxz(1)); // frdp->fdp
// 		dp.ifftxz();                  // fdp->dp
// 	gettimeofday(&time1, NULL); t_getdp += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);

// 	gettimeofday(&time0, NULL);
// 		// upcalc
// 		#pragma omp parallel
// 		update(fld, fldh, dt);
// 	gettimeofday(&time1, NULL); t_update += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
		
// 	if (cnt == 10) {
// 		cout << endl << endl
// 			 << "total:\t"    << 1./cnt * (t_urhs+t_getuh+t_getdp+t_update) << endl
// 			 << "t_getuh:\t"  << t_getuh / t_urhs << endl
// 			 << "t_getdp:\t"  << t_getdp / t_urhs << endl
// 			 << "t_update:\t" << t_update/ t_urhs << endl << endl;
// 		exit(0);
// 	}
// }




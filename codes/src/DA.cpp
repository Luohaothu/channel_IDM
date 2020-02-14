# include <iostream>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <cmath>

# include "DA.h"
# include "Matrix.h"

using namespace std;




bool DA::getExp(double time, const Vctr &UE)
{
	if (_iter) {
		setMesh(_FLDH.meshGet());
		for (int j=1; j<Ny; j++) {
		for (int k=0; k<Nz; k++) {
		for (int i=0; i<Nx; i++) {
			_FLDH.S.id(i,j,k) = _F.module(i,j,k);
		}}}
		cout << "Time: "     << time
		     << "\tIter: "   << _iter
		     << "\tError: "  << _erro
		     << "\tMean F: " << _FLDH.S.bulkMeanU() << endl;
	}

	if (_iter) _F.reset(_iter = 0);

	double interval = .5;
	if (fabs(fmod(time+INFTSM/2., interval)) < INFTSM) {
		_UE[1] = UE[1];
		_UE[2] = UE[2];
		_UE[3] = UE[3];
		return true;
	}
	return false; // no experiment data exist
}

bool DA::ifIter(const Vctr &U, double e, int n)
{
	if (_iter >= n) return false;
	
	((_FLDH.V[1] = U[1]) -= _UE[1]) *= _MSK[1];
	((_FLDH.V[2] = U[2]) -= _UE[2]) *= _MSK[2];
	((_FLDH.V[3] = U[3]) -= _UE[3]) *= _MSK[3];

	setMesh(U.meshGet());
	for (int j=1; j<Ny; j++) {
	for (int k=0; k<Nz; k++) {
	for (int i=0; i<Nx; i++) {
		_FLDH.S.id(i,j,k) = _FLDH.V.module(i,j,k) / _UE.module(i,j,k);
	}}}
	return (_erro = _FLDH.S.bulkMeanU() / _vexp) >= e;
}

const Vctr& DA::getForce(double alpha)
{
	double lambda = 0.;

	setMesh(_FLDH.meshGet());
	for (int j=1; j<Ny; j++) {
	for (int k=0; k<Nz; k++) {
	for (int i=0; i<Nx; i++) {
		lambda = fmax(lambda, _FLDH.V.module(i,j,k));
	}}}

	lambda = alpha / lambda;

	_F[1] -= (_FLDH.V[1] *= lambda); // NOTE: in Liu's paper, the sign before V should be +
	_F[2] -= (_FLDH.V[2] *= lambda);
	_F[3] -= (_FLDH.V[3] *= lambda);
	_iter ++;

	return _F;
}

void DA::getAdj(const Vctr &U, const Feld &VIS, double dt)
{
	urhs(_FLDH.V, U, _UE, _MSK);

	getuh1(_FLDH.V, U, VIS, dt);
	getuh2(_FLDH.V, U, VIS, dt);
	getuh3(_FLDH.V, U, VIS, dt);

	rhsdp(_FLDH.S, _FLDH.V, dt);
	_FLDH.S.fft();
	getfdp(_FLDH.S, 0.);
	_FLDH.S.ifft();

	update(_FLDH.V, _FLDH.S, dt);

	// homogeneous BC applied automatically
}


/***** adjoint velocity computation *****/

void DA::urhs(Vctr &UH, const Vctr &U, const Vctr &UE, const Vctr &MSK)
/* compute the RHS (without mbc) of momentum equation and store in UH */
{
	(((UH[1] = U[1]) -= UE[1]) *= MSK[1]) *= 2.;
	(((UH[2] = U[2]) -= UE[2]) *= MSK[2]) *= 2.;
	(((UH[3] = U[3]) -= UE[3]) *= MSK[3]) *= 2.;

	// mbc, cbc = 0 for homogeneous BC
}

void DA::getuh1(Vctr &UH, const Vctr &U, const Feld &VIS, double dt)
/* compute deltaU^**, result returned by uh (the RHS ruh should be pre stored in uh ) */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, imjp, imkp;
	double vis1, vis2, vis3, vis4, vis5, vis6, a;

	const Mesh &ms = setMesh(UH.meshGet());
	double *u,*v,*w;             U.ptrGet(u,v,w);
	double *uh,*vh,*wh;         UH.ptrGet(uh,vh,wh);
	double *nux,*nuy,*nuz,*nu; VIS.ptrGet(nux,nuy,nuz,nu);

	double *api = new double [Nx], *aci = new double [Nx], *ami = new double [Nx], *R1 = new double [Nx];
	double *apj = new double [Ny], *acj = new double [Ny], *amj = new double [Ny], *R2 = new double [Ny];
	double *apk = new double [Nz], *ack = new double [Nz], *amk = new double [Nz], *R3 = new double [Nz];
	Matrix matx(Nx);
	Matrix maty(Ny-1);
	Matrix matz(Nz);

	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		for (j=1; j<Ny; j++) {
			idx = ms.IDX(i,j,k);
			ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
			im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);
			imjp = ms.IDX(ima[i],j+1,k);

			a = .25 * (v[idx] + v[im] + v[jp] + v[imjp]);
			vis3 = nuz[idx];
			vis4 = nuz[jp];

			apj[j] = (-.5/h[j+1]*a                                  - 1./dy[j] * vis4/h[j+1]            ) * dt;
			acj[j] = (-.5/dy[j] *a *(dy[j+1]/h[j+1] - dy[j-1]/h[j]) + 1./dy[j] *(vis4/h[j+1]+vis3/h[j]) ) * dt + 1;
			amj[j] = ( .5/h[j]  *a                                  - 1./dy[j] * vis3/h[j]              ) * dt;

			R2 [j] = dt * uh[idx];
		}
		maty.tdma( & amj[1], & acj[1], & apj[1], & R2[1] ); // apj at j=Ny-1 and amj at j=1 are redundant in tdma, thus no need for explicit removal
		for (j=1; j<Ny; j++) uh[ms.IDX(i,j,k)] = R2[j];
	}}

	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			idx = ms.IDX(i,j,k);
			ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
			im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);

			a = .5/dx * (u[ip] - u[im]);
			vis1 = nu[im];
			vis2 = nu[idx];

			api[i] = ( .5/dx * u[idx] - 2./dx2 * vis2       ) * dt;
			aci[i] = ( a              + 2./dx2 *(vis2+vis1) ) * dt + 1;
			ami[i] = (-.5/dx * u[idx] - 2./dx2 * vis1       ) * dt;

			R1 [i] = uh[idx];
		}
		matx.ctdma( ami, aci, api, R1 );
		for (i=0; i<Nx; i++) uh[ms.IDX(i,j,k)] = R1[i];
	}}

	for (j=1; j<Ny; j++) {
	for (i=0; i<Nx; i++) {
		for (k=0; k<Nz; k++) {
			idx = ms.IDX(i,j,k);
			ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
			im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);
			imkp = ms.IDX(ima[i],j,kpa[k]);

			a = .25 * (w[idx] + w[im] + w[kp] + w[imkp]);
			vis5 = nuy[idx];
			vis6 = nuy[kp];

			apk[k] = (-.5/dz * a - 1./dz2 * vis6        ) * dt;
			ack[k] = (             1./dz2 *(vis6+vis5)  ) * dt + 1;
			amk[k] = ( .5/dz * a - 1./dz2 * vis5        ) * dt;

			R3 [k] = uh[idx];
		}
		matz.ctdma( amk, ack, apk, R3 );
		for (k=0; k<Nz; k++) uh[ms.IDX(i,j,k)] = R3[k];
	}}

	delete [] api; delete [] aci; delete [] ami; delete [] R1;
	delete [] apj; delete [] acj; delete [] amj; delete [] R2;
	delete [] apk; delete [] ack; delete [] amk; delete [] R3;
}


void DA::getuh2(Vctr &UH, const Vctr &U, const Feld &VIS, double dt)
/* compute deltaU^**, result returned by vh (the RHS rvh should be pre stored in uh ) */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, ipjm, jmkp, imjm, jmkm;
	double vis1, vis2, vis3, vis4, vis5, vis6, u1, u2, a, m21uh;

	const Mesh &ms = setMesh(UH.meshGet());
	double *u,*v,*w;             U.ptrGet(u,v,w);
	double *uh,*vh,*wh;         UH.ptrGet(uh,vh,wh);
	double *nux,*nuy,*nuz,*nu; VIS.ptrGet(nux,nuy,nuz,nu);

	double *api = new double [Nx], *aci = new double [Nx], *ami = new double [Nx], *R1 = new double [Nx];
	double *apj = new double [Ny], *acj = new double [Ny], *amj = new double [Ny], *R2 = new double [Ny];
	double *apk = new double [Nz], *ack = new double [Nz], *amk = new double [Nz], *R3 = new double [Nz];
	Matrix matx(Nx);
	Matrix maty(Ny-2);
	Matrix matz(Nz);

	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		for (j=2; j<Ny; j++) {
			idx = ms.IDX(i,j,k);
			ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
			im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);
			ipjm = ms.IDX(ipa[i],j-1,k); imjm = ms.IDX(ima[i],j-1,k);

			a = .5/h[j] * (v[jp] - v[jm]);
			vis3 = nu[jm];
			vis4 = nu[idx];

			apj[j] = (-.5/h[j] * v[idx] - 2./h[j] * vis4/dy[j]               ) * dt;
			acj[j] = ( a                + 2./h[j] *(vis4/dy[j]+vis3/dy[j-1]) ) * dt + 1;
			amj[j] = ( .5/h[j] * v[idx] - 2./h[j] * vis3/dy[j-1]             ) * dt;

			// m21uh
			a = .5/dx * (v[ip] - v[im]);
			vis1 = nuz[idx];
			vis2 = nuz[ip];
			u2 = .5/h[j] * (uh[ip]*dy[j-1] + uh[ipjm]*dy[j]);
			u1 = .5/h[j] * (uh[idx]*dy[j-1] + uh[jm]*dy[j]);

			m21uh = .5*a * (u1+u2) - 1./dx/h[j] * ( vis2 * (uh[ip]-uh[ipjm]) - vis1 * (uh[idx]-uh[jm]) );

			R2 [j] = dt * ( vh[idx] - m21uh );
		}
		maty.tdma( & amj[2], & acj[2], & apj[2], & R2[2] ); // apj at j=Ny-1 and amj at j=1 are redundant in tdma, thus no need for explicit removal
		for (j=2; j<Ny; j++) vh[ms.IDX(i,j,k)] = R2[j];
	}}

	for (j=2; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			idx = ms.IDX(i,j,k);
			ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
			im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);
			ipjm = ms.IDX(ipa[i],j-1,k); imjm = ms.IDX(ima[i],j-1,k);

			a = .25/h[j] * ( (u[idx]+u[ip]) * dy[j-1] + (u[jm]+u[ipjm]) * dy[j] );
			vis1 = nuz[idx];
			vis2 = nuz[ip];

			api[i] = (-.5/dx * a - 1./dx2 * vis2       ) * dt;
			aci[i] = (             1./dx2 *(vis2+vis1) ) * dt + 1;
			ami[i] = ( .5/dx * a - 1./dx2 * vis1       ) * dt;

			R1 [i] = vh[idx];
		}
		matx.ctdma( ami, aci, api, R1 );
		for (i=0; i<Nx; i++) vh[ms.IDX(i,j,k)] = R1[i];
	}}

	for (j=2; j<Ny; j++) {
	for (i=0; i<Nx; i++) {
		for (k=0; k<Nz; k++) {
			idx = ms.IDX(i,j,k);
			ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
			im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);
			jmkp = ms.IDX(i,j-1,kpa[k]); jmkm = ms.IDX(i,j-1,kma[k]);

			a = .25/h[j] * ( (w[idx]+w[kp]) * dy[j-1] + (w[jm]+w[jmkp]) * dy[j] );
			vis5 = nux[idx];
			vis6 = nux[kp];

			apk[k] = (-.5/dz * a - 1./dz2 * vis6       ) * dt;
			ack[k] = (             1./dz2 *(vis6+vis5) ) * dt + 1;
			amk[k] = ( .5/dz * a - 1./dz2 * vis5       ) * dt;

			R3 [k] = vh[idx];
		}
		matz.ctdma( amk, ack, apk, R3 );
		for (k=0; k<Nz; k++) vh[ms.IDX(i,j,k)] = R3[k];
	}}
	
	delete [] api; delete [] aci; delete [] ami; delete [] R1;
	delete [] apj; delete [] acj; delete [] amj; delete [] R2;
	delete [] apk; delete [] ack; delete [] amk; delete [] R3;
}


void DA::getuh3(Vctr &UH, const Vctr &U, const Feld &VIS, double dt)
/* compute deltaW^*, deltaV^*, deltaU^*, and update U^*, V^*, W^* (the RHS rwh should be pre stored in wh ) */
{
	int i, j, k, idx, ip, im, jp, jm, kp, km, jup, jum;
	int ipkm, jpkm, imjp, imkp, jmkp, imkm, jmkm, imjm;
	double vis1, vis2, vis3, vis4, vis5, vis6;
	double u1, u2, v1, v2, w1, w2, a;
	double m31uh, m32vh, m23wh, m12vh, m13wh;

	const Mesh &ms = setMesh(UH.meshGet());
	double *u,*v,*w;             U.ptrGet(u,v,w);
	double *uh,*vh,*wh;         UH.ptrGet(uh,vh,wh);
	double *nux,*nuy,*nuz,*nu; VIS.ptrGet(nux,nuy,nuz,nu);

	double *api = new double [Nx], *aci = new double [Nx], *ami = new double [Nx], *R1 = new double [Nx];
	double *apj = new double [Ny], *acj = new double [Ny], *amj = new double [Ny], *R2 = new double [Ny];
	double *apk = new double [Nz], *ack = new double [Nz], *amk = new double [Nz], *R3 = new double [Nz];
	Matrix matx(Nx);
	Matrix maty(Ny-1);
	Matrix matz(Nz);

	// ( I + dt M_33^2 )
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		for (j=1; j<Ny; j++) {  jup = (j!=Ny-1); jum = (j!=1);
			idx = ms.IDX(i,j,k);
			ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
			im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);
			ipkm = ms.IDX(ipa[i],j,kma[k]); jpkm = ms.IDX(i,j+1,kma[k]);
			imkm = ms.IDX(ima[i],j,kma[k]); jmkm = ms.IDX(i,j-1,kma[k]);

			a = .25 * (v[idx] + v[jp] + v[km] + v[jpkm]);
			vis3 = nux[idx];
			vis4 = nux[jp];

			apj[j] = (-.5/h[j+1]* a                                  - 1./dy[j] * vis4/h[j+1]            ) * dt;
			acj[j] = (-.5/dy[j] * a *(dy[j+1]/h[j+1] - dy[j-1]/h[j]) + 1./dy[j] *(vis4/h[j+1]+vis3/h[j]) ) * dt + 1;
			amj[j] = ( .5/h[j]  * a                                  - 1./dy[j] * vis3/h[j]              ) * dt;

			// m31uh
			a = .5/dx * (w[ip] - w[im]);
			vis1 = nuy[idx];
			vis2 = nuy[ip];
			u2 = .5 * (uh[ip] + uh[ipkm]);
			u1 = .5 * (uh[idx] + uh[km]);

			m31uh = .5*a * (u2+u1) + .1/dx/dz * ( vis2 * (uh[ip]-uh[ipkm]) - vis1 * (uh[idx]-uh[km]) );

			// m32vh
			w2 = .5/h[j+1] * (w[idx]*dy[j+1] + w[jp]*dy[j]); //* jup;
			w1 = .5/h[j]   * (w[idx]*dy[j-1] + w[jm]*dy[j]); //* jum;
			a = 1./dy[j] * (w2-w1);
			vis3 = nux[idx] * jum;
			vis4 = nux[jp]  * jup;
			v2 = .5 * (vh[jp] + vh[jpkm]);
			v1 = .5 * (vh[idx] + vh[km]);

			m32vh = .5*a * (v2+v1) - 1./dy[j]/dz * ( vis4 * (vh[jp]-vh[jpkm]) - vis3 * (vh[idx]-vh[km]) );

			R2 [j] = dt * ( wh[idx] - m31uh - m32vh );
		}
		maty.tdma( & amj[1], & acj[1], & apj[1], & R2[1] );
		for (j=1; j<Ny; j++) wh[ms.IDX(i,j,k)] = R2[j];
	}}

	// ( I + dt M_33^1 )
	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
		for (i=0; i<Nx; i++) {
			idx = ms.IDX(i,j,k);
			ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
			im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);
			ipkm = ms.IDX(ipa[i],j,kma[k]); imkm = ms.IDX(ima[i],j,kma[k]);

			a = .25 * (u[idx] + u[ip] + u[km] + u[ipkm]);
			vis1 = nuy[idx];
			vis2 = nuy[ip];

			api[i] = (-.5/dx * a - 1./dx2 * vis2       ) * dt;
			aci[i] = (             1./dx2 *(vis2+vis1) ) * dt + 1;
			ami[i] = ( .5/dx * a - 1./dx2 * vis1       ) * dt;

			R1 [i] = wh[idx];
		}
		matx.ctdma( ami, aci, api, R1 );
		for (i=0; i<Nx; i++) wh[ms.IDX(i,j,k)] = R1[i];
	}}

	// ( I + dt M_33^3 )
	for (j=1; j<Ny; j++) {
	for (i=0; i<Nx; i++) {
		for (k=0; k<Nz; k++) {
			idx = ms.IDX(i,j,k);
			ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
			im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);

			a = .5/dz * (w[kp]-w[km]);
			vis5 = nu[km];
			vis6 = nu[idx];

			apk[k] = (-.5/dz * w[idx] - 2./dz2 * vis6       ) * dt;
			ack[k] = ( a              + 2./dz2 *(vis6+vis5) ) * dt + 1;
			amk[k] = ( .5/dz * w[idx] - 2./dz2 * vis5       ) * dt;

			R3 [k] = wh[idx];
		}
		matz.ctdma( amk, ack, apk, R3 );
		for (k=0; k<Nz; k++) wh[ms.IDX(i,j,k)] = R3[k];
	}}

	delete [] api; delete [] aci; delete [] ami; delete [] R1;
	delete [] apj; delete [] acj; delete [] amj; delete [] R2;
	delete [] apk; delete [] ack; delete [] amk; delete [] R3;

	// update dvh
	for (j=2; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = ms.IDX(i,j,k);
		ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
		im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);
		jmkp = ms.IDX(i,j-1,kpa[k]); jmkm = ms.IDX(i,j-1,kma[k]);

		// m23wh
		a = .5/dz * (v[kp]-v[km]);
		w2 = .5/h[j] * (wh[kp]*dy[j-1] + wh[jmkp]*dy[j]);
		w1 = .5/h[j] * (wh[idx]*dy[j-1] + wh[jm]*dy[j]);
		vis5 = nux[idx];
		vis6 = nux[kp];

		m23wh = .5*a * (w2+w1) - 1./dz/h[j] * ( vis6 * (wh[kp]-wh[jmkp]) - vis5 * (wh[idx]-wh[jm]) );

		vh[idx] -= dt * m23wh;
	}}}
	
	// update duh
	for (j=1; j<Ny; j++) {  jup = (j!=Ny-1); jum = (j!=1);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = ms.IDX(i,j,k);
		ip = ms.IDX(ipa[i],j,k); jp = ms.IDX(i,j+1,k); kp = ms.IDX(i,j,kpa[k]);
		im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);
		imjp = ms.IDX(ima[i],j+1,k); imkp = ms.IDX(ima[i],j,kpa[k]);
		imjm = ms.IDX(ima[i],j-1,k); imkm = ms.IDX(ima[i],j,kma[k]);

		// m12vh
		u2 = .5/h[j+1] * (u[idx]*dy[j+1] + u[jp]*dy[j]); //* jup;
		u1 = .5/h[j]   * (u[idx]*dy[j-1] + u[jm]*dy[j]); //* jum;
		a = 1./dy[j] * (u2-u1);
		vis3 = nuz[idx] * jum;
		vis4 = nuz[jp]  * jup;
		v2 = .5 * (vh[jp] + vh[imjp]);
		v1 = .5 * (vh[idx] + vh[im]);

		m12vh = .5*a * (v2+v1) - 1./dy[j]/dx * ( vis4 * (vh[jp]-vh[imjp]) - vis3 * (vh[idx]-vh[im]) );

		// m13wh
		a = .5/dz * (u[kp]-u[km]);
		w2 = .5 * (wh[kp] + wh[imkp]);
		w1 = .5 * (wh[idx] + wh[im]);
		vis5 = nuy[idx];
		vis6 = nuy[kp];

		m13wh = .5*a * (w2+w1) - 1./dz/dx * ( vis6 * (wh[kp]-wh[imkp]) - vis5 * (wh[idx]-wh[im]) );

		uh[idx] -= dt * (m12vh + m13wh);
	}}}
}



/***** projector computation *****/

void DA::rhsdp(Scla &DP, const Vctr &UH, double dt)
/* compute RHS of Poisson equation and store in dp */
{
	int i, j, k, idx, ip, jp, kp, j1, j0, jup, jum;

	const Mesh &ms = setMesh(DP.meshGet());
	double *uh,*vh,*wh; UH.ptrGet(uh,vh,wh);

	for (j=1; j<Ny; j++) { jup = (j!=Ny-1); jum = (j!=1);
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx= ms.IDX(i,j,k); ip = ms.IDX(ipa[i],j,k);
		j1 = ms.IDX(i,1,k); jp = ms.IDX(i,j+1,k);
		j0 = ms.IDX(i,0,k); kp = ms.IDX(i,j,kpa[k]);

		// ( Du^h - cbc ) / dt
		DP.id(idx) = 1.0/dt * (
			( uh[ip] - uh[idx] ) / dx
		+	( wh[kp] - wh[idx] ) / dz
		+	( vh[jp] * jup //+ vbc[j1] * (1-jup) // mbc, cbc = 0 for homogeneous BC
			- vh[idx]* jum //- vbc[j0] * (1-jum)
			                   ) / dy[j]	);
	}}}
}

void DA::getfdp(Scla &DP, double refp)
/* compute FDP (in Fourier space), result returned by fdp (the RHS frdp should be pre stored in fdp ) */
{
	int i, j, k, Nxc = (int) (Nx/2+1);

	const Mesh &ms = setMesh(DP.meshGet());
	double *cpj = new double [Ny];
	double *ccj = new double [Ny];
	double *cmj = new double [Ny];
	double *cfj1= new double [Ny], *cfj2= new double [Ny];
	Matrix matyp(Ny-1);

	for (k=0; k<Nz; k++) {
	for (i=0; i<Nxc; i++) { // negative k_x need not solve
		for (j=1; j<Ny; j++) {
			cpj[j] = ms.ppj[j];
			ccj[j] = ms.pcj[j] - ms.ak3[k] - ms.ak1[i];
			cmj[j] = ms.pmj[j];
			cfj1[j] = DP.idf(2*i,  j,k); // real part
			cfj2[j] = DP.idf(2*i+1,j,k); // imaginary part
		}
		// set reference pressure P(kx=0,kz=0,j=1) to 0
		if (k==0 && i==0) {
			cpj[1] = 0;
			ccj[1] = 1;
			cmj[1] = 0;
			cfj1[1] = - Nxz * refp;	// fdp(kx=0,kz=0,j=1) = -Nxz * fp(kx=0,kz=0,j=1) ==> <dp(j=1)> = - <p(j=1)>
			cfj2[1] = 0;
		}

		matyp.tdma( & cmj[1], & ccj[1], & cpj[1], & cfj1[1] );
		matyp.tdma( & cmj[1], & ccj[1], & cpj[1], & cfj2[1] );

		for (j=1; j<Ny; j++) {
			DP.idf(2*i,  j,k) = cfj1[j];
			DP.idf(2*i+1,j,k) = cfj2[j];
		}
	}}

	delete [] cpj; delete [] ccj; delete [] cmj;
	delete [] cfj1;delete [] cfj2;
}

void DA::update(Vctr &UH, const Scla &DP, double dt)
/* project from U^* to U^n+1 using DP */
{
	int i, j, k, idx, im, jm, km;
	const Mesh &ms = setMesh(UH.meshGet());
	double *uh,*vh,*wh; UH.ptrGet(uh,vh,wh);
	double *dp = DP.blkGet();

	for (j=1; j<Ny; j++) {
	for (k=0; k<Nz; k++) {
	for (i=0; i<Nx; i++) {
		idx = ms.IDX(i,j,k);
		im = ms.IDX(ima[i],j,k); jm = ms.IDX(i,j-1,k); km = ms.IDX(i,j,kma[k]);

		uh[idx] -= dt/dx * (dp[idx] - dp[im]); if (j>1)
		vh[idx] -= dt/h[j]*(dp[idx] - dp[jm]);
		wh[idx] -= dt/dz * (dp[idx] - dp[km]);
	}}}
}




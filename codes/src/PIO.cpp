#include "PIO.h"
#include "Interp.h"
#include "Filter.h"
#include "Bcond.h"


using namespace std;

// #define MODULATION


void read_PIO_Sup(double *hr, double *hi,
	double &alfv, double &dxav, double &dzav,
	double &alfw, double &dxaw, double &dzaw,
	double lc, double y, const double *kx, int Nxc);

void read_PIO_Mod(
	double &gmaup, double &dxup, double &dzup,
	double &gmaum, double &dxum, double &dzum,
	double &gmav,  double &dxv,  double &dzv,
	double &gmaw,  double &dxw,  double &dzw,
	double lc, double y);


void PIO::PredictBoundary(Vctr &veltmp, const Vctr &vel, const Vctr &velmfu, double Ret, double rsclx, double rsclu)
{
	const Mesh &ms = vel.ms;
	const Mesh &ms0 = velmfu.ms;

	double yb1 = ms.y(1); // position to predict
	double yb2 = ms.y(ms.Ny);

	double yo1 = 3.9 / sqrt(Ret); // position of outer signal
	double yo2 = 2. - yo1;

	Vctr &velo = veltmp;
	Vctr &velb = veltmp;
	Scla &uo = velo[1], &vo = velo[2], &wo = velo[3]; // to hold uL
	Scla &ub = velb[1], &vb = velb[2], &wb = velb[3]; // to hold u_p

	// PIO coefficients
	double *hr = new double[ms.Nxc], alfv, dxav, dzav;
	double *hi = new double[ms.Nxc], alfw, dxaw, dzaw;

	read_PIO_Sup(hr, hi, alfv, dxav, dzav,
					 alfw, dxaw, dzaw, 1./Ret, yb1, ms.kx(), ms.Nxc);

	double gmaup, dxup, dzup;
	double gmaum, dxum, dzum;
	double gmav,  dxv,  dzv;
	double gmaw,  dxw,  dzw;

	read_PIO_Mod(gmaup, dxup, dzup,
				 gmaum, dxum, dzum,
				 gmav,  dxv,  dzv,
				 gmaw,  dxw,  dzw, 1./Ret, yb1);


	int jo1 = 2; // position in veltmp to store uOL
	int jo2 = jo1 + 1;
	
	#pragma omp parallel
	{
		// interpolate from LES outer region, mean values will be removed later
		#pragma omp for collapse(2)
		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			double x = ms.x(i), xc = ms.xc(i);
			double z = ms.z(k), zc = ms.zc(k);

			uo(i,jo1,k) = Interp::InterpNodeU(x, yo1,zc,vel[1]);
			vo(i,jo1,k) = Interp::InterpNodeV(xc,yo1,zc,vel[2]);
			wo(i,jo1,k) = Interp::InterpNodeW(xc,yo1,z, vel[3]);

			uo(i,jo2,k) = Interp::InterpNodeU(x, yo2,zc,vel[1]);
			vo(i,jo2,k) = Interp::InterpNodeV(xc,yo2,zc,vel[2]);
			wo(i,jo2,k) = Interp::InterpNodeW(xc,yo2,z, vel[3]);
		}}

		uo.fftx (jo1, jo2);
		vo.fftxz(jo1, jo2);
		wo.fftxz(jo1, jo2);

	// // keep the random perturbation of vo constant in time and random in space
	// static bool randinput = true;
	// static double randtable[1024*1024];
	// if (randinput) {
	// 	FILE* fp = fopen("randtable", "rb");
	// 	fread(randtable, sizeof(double), 1024*1024, fp);
	// 	fclose(fp);
	// 	randinput = false;
	// }

		#pragma omp for collapse(2)
		for (int j=jo1; j<=jo2; j++) {
		for (int i=0; i<ms.Nxc; i++) {

			// get mean value at wave number k_x = 0
			double um = 0;
			if (i == 0)
				for (int k=1; k<ms.Nz; k++)
					um += uo(2*i,j,k) / (ms.Nz-1.);

			for (int k=1; k<ms.Nz; k++) {
				double ur = uo(2*i,  j,k) - um;
				double ui = uo(2*i+1,j,k);
				// H needs conj because the fft defined in calibration is inversed
				uo(2*i,  j,k) = hr[i] * ur + hi[i] * ui;
				uo(2*i+1,j,k) = hr[i] * ui - hi[i] * ur;
			}
		}}

		#pragma omp for collapse(2)
		for (int k=0; k<ms.Nz-1; k++) {
		for (int i=0; i<ms.Nxc; i++) {
			// keep only the fluctuations with scale strictly larger than MFU
			if (fabs(ms.kx(i)/rsclx) >= 2.*PI/ms0.Lx ||
				fabs(ms.kz(k)/rsclx) >= 2.*PI/ms0.Lz ||
				(k==0 && i==0))
			{
				vo(2*i,jo1,k) = 0; vo(2*i+1,jo1,k) = 0;
				wo(2*i,jo1,k) = 0; wo(2*i+1,jo1,k) = 0;
				vo(2*i,jo2,k) = 0; vo(2*i+1,jo2,k) = 0;
				wo(2*i,jo2,k) = 0; wo(2*i+1,jo2,k) = 0;
			}
			// // offset vOL by half channel width (and/or length) to break the correlation between u & v
			// else if (k%2 ^ i%2) {
			// 	vo(2*i,1,k) *= -1; vo(2*i+1,1,k) *= -1;
			// 	vo(2*i,2,k) *= -1; vo(2*i+1,2,k) *= -1;
			// } // test shows that <uLvL> is not sufficiently eliminated in this way
			// // randomly perturb the phase (keeping energy spectra unchanged) to interrupt the correlation between u & v
			// else {
			// 	double phi=k*PI, vr, vi;
				
			// 	phi = 2*PI * randtable[k*mso.Nxc+i];
			// 	vr = vo(2*i,  1,k);
			// 	vi = vo(2*i+1,1,k);

			// 	vo(2*i,  1,k) = vr * cos(phi) - vi * sin(phi);
			// 	vo(2*i+1,1,k) = vr * sin(phi) + vi * cos(phi);

			// 	phi = 2*PI * randtable[k*mso.Nxc+i+(mso.Nz-1)*mso.Nxc];
			// 	vr = vo(2*i,  2,k);
			// 	vi = vo(2*i+1,2,k);

			// 	vo(2*i,  2,k) = vr * cos(phi) - vi * sin(phi);
			// 	vo(2*i+1,2,k) = vr * sin(phi) + vi * cos(phi);

			// 	// if (i==10 && k==10) cout << "1 " << phi << endl;
			// 	// if (i==5 && k==5) cout << "0 " << phi << endl;
			// }
		}}

		uo.ifftx (jo1, jo2);
		vo.ifftxz(jo1, jo2);
		wo.ifftxz(jo1, jo2);
	}

	Bcond::SetBoundaryX(velo);
	Bcond::SetBoundaryZ(velo);

	delete[] hr;
	delete[] hi;


	// // check large-scale Reynolds stress
	// double uv = 0;
	// for (int k=1; k<mso.Nz; k++) {
	// for (int i=1; i<mso.Nx; i++) {
	// 	uv += uo(i,0,k) * vo(i,1,k) - uo(i,2,k) * vo(i,2,k);
	// }}
	// uv /= 2. * (mso.Nx-1) * (mso.Nz-1);
	// cout << uv << endl;


	// filter from MFU matching inner scale
	for (int j=0; j<=ms.Ny; j+=(j ? ms.Ny-1 : 1)) {
	#pragma omp parallel for collapse(2)
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		double x  = ms.x (i) * rsclx, xc = ms.xc(i) * rsclx;
		double z  = ms.z (k) * rsclx, zc = ms.zc(k) * rsclx;
		double dx = ms.dx(i) * rsclx, hx = ms.hx(i) * rsclx;
		double dz = ms.dz(k) * rsclx, hz = ms.hz(k) * rsclx;
		double y  = Filter::WallRscl(ms.y (j), rsclx);
		double yc = Filter::WallRscl(ms.yc(j), rsclx);

		if (i>0) ub(i,j,k) = Filter::FilterNodeU(x,yc,zc, hx,0,dz, velmfu[1]) * rsclu;
		if (j>0) vb(i,j,k) = Filter::FilterNodeV(xc,y,zc, dx,0,dz, velmfu[2]) * rsclu;
		if (k>0) wb(i,j,k) = Filter::FilterNodeW(xc,yc,z, dx,0,hz, velmfu[3]) * rsclu;

		// handle HALFMFU
		if (1+ms.Ny < 2/ms0.Ly*j) vb(i,j,k) *= -1;
	}}}

#ifdef MODULATION
	// modulation
	for (int j=0; j<=ms.Ny; j+=(j ? ms.Ny-1 : 1)) {

		double um = ub.meanxz(j);
		double vm = vb.meanxz(j);
		double wm = wb.meanxz(j);
		double uom= uo.meanxz(j<=1 ? jo1 : jo2);

		#pragma omp parallel for collapse(2)
		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			double x = ms.x(i), xc = ms.xc(i);
			double z = ms.z(k), zc = ms.zc(k);
			double y = j<=1 ? ms.yc(jo1) : ms.yc(jo2);

			double modu;
			double modv = gmav * (Interp::InterpNodeU(xc+dxv,y,zc+dzv,uo) - uom);
			double modw = gmaw * (Interp::InterpNodeU(xc+dxw,y,z +dzw,uo) - uom);

			if (ub(i,j,k) > um) modu = gmaup * (Interp::InterpNodeU(x+dxup,y,zc+dzup,uo) - uom);
			else                modu = gmaum * (Interp::InterpNodeU(x+dxum,y,zc+dzum,uo) - uom);

			if (j != 1) ub(i,j,k) = um + (ub(i,j,k) - um) * (1 + modu);
			if (j != 0) vb(i,j,k) = vm + (vb(i,j,k) - vm) * (1 + modv);
			if (j != 1) wb(i,j,k) = wm + (wb(i,j,k) - wm) * (1 + modw);
		}}
	}
#endif

	// combine small- & large-scales
	#pragma omp parallel for collapse(2)
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		ub(i,0,    k) += uo(i,jo1,k);
		ub(i,ms.Ny,k) += uo(i,jo2,k);

		double x = ms.xc(i) + dxav;
		double z = ms.zc(k) + dzav;

		vb(i,1,    k) += alfv * Interp::InterpNodeV(x,ms.y(jo1),z,vo);
		vb(i,ms.Ny,k) += alfv * Interp::InterpNodeV(x,ms.y(jo2),z,vo);

		x = ms.xc(i) + dxaw;
		z = ms.z (k) + dzaw;

		wb(i,0,    k) += alfw * Interp::InterpNodeW(x,ms.yc(jo1),z,wo);
		wb(i,ms.Ny,k) += alfw * Interp::InterpNodeW(x,ms.yc(jo2),z,wo);
	}}

	Bcond::SetBoundaryX(velb);
	Bcond::SetBoundaryZ(velb);

	// record for debug
	// static int cnt = 0;
	// if ((++cnt) % 5000 == 0)
	// 	PIO::Predict(yb1, vel, velmfu, Ret, rsclx, rsclu, cnt);
}

void PIO::PredictBoundary(Vctr &veltmp, const Vctr &vel, double time, double Ret)
{
	const Mesh &ms = vel.ms;

	double yb1 = ms.y(1); // position to predict
	double yb2 = ms.y(ms.Ny);

	double yo1 = 3.9 / sqrt(Ret); // position of outer signal
	double yo2 = 2. - yo1;

	Vctr &velo = veltmp;
	Vctr &velb = veltmp;
	Scla &uo = velo[1], &vo = velo[2], &wo = velo[3]; // to hold uL
	Scla &ub = velb[1], &vb = velb[2], &wb = velb[3]; // to hold u_p

	// PIO coefficients
	double *hr = new double[ms.Nxc], alfv, dxav, dzav;
	double *hi = new double[ms.Nxc], alfw, dxaw, dzaw;

	read_PIO_Sup(hr, hi, alfv, dxav, dzav,
					 alfw, dxaw, dzaw, 1./Ret, yb1, ms.kx(), ms.Nxc);

	double gmaup, dxup, dzup;
	double gmaum, dxum, dzum;
	double gmav,  dxv,  dzv;
	double gmaw,  dxw,  dzw;

	read_PIO_Mod(gmaup, dxup, dzup,
				 gmaum, dxum, dzum,
				 gmav,  dxv,  dzv,
				 gmaw,  dxw,  dzw, 1./Ret, yb1);


	int jo1 = 2; // position in veltmp to store uOL
	int jo2 = jo1 + 1;
	
	#pragma omp parallel
	{
		// interpolate from LES outer region, mean values will be removed later
		#pragma omp for collapse(2)
		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			double x = ms.x(i), xc = ms.xc(i);
			double z = ms.z(k), zc = ms.zc(k);

			uo(i,jo1,k) = Interp::InterpNodeU(x, yo1,zc,vel[1]);
			vo(i,jo1,k) = Interp::InterpNodeV(xc,yo1,zc,vel[2]);
			wo(i,jo1,k) = Interp::InterpNodeW(xc,yo1,z, vel[3]);

			uo(i,jo2,k) = Interp::InterpNodeU(x, yo2,zc,vel[1]);
			vo(i,jo2,k) = Interp::InterpNodeV(xc,yo2,zc,vel[2]);
			wo(i,jo2,k) = Interp::InterpNodeW(xc,yo2,z, vel[3]);
		}}

		uo.fftx (jo1, jo2);
		vo.fftxz(jo1, jo2);
		wo.fftxz(jo1, jo2);

		#pragma omp for collapse(2)
		for (int j=jo1; j<=jo2; j++) {
		for (int i=0; i<ms.Nxc; i++) {

			// get mean value at wave number k_x = 0
			double um = 0;
			if (i == 0)
				for (int k=1; k<ms.Nz; k++)
					um += uo(2*i,j,k) / (ms.Nz-1.);

			for (int k=1; k<ms.Nz; k++) {
				double ur = uo(2*i,  j,k) - um;
				double ui = uo(2*i+1,j,k);
				// H needs conj because the fft defined in calibration is inversed
				uo(2*i,  j,k) = hr[i] * ur + hi[i] * ui;
				uo(2*i+1,j,k) = hr[i] * ui - hi[i] * ur;
			}
		}}

		#pragma omp for collapse(2)
		for (int k=0; k<ms.Nz-1; k++) {
		for (int i=0; i<ms.Nxc; i++) {
			// keep only the fluctuations with scale strictly larger than MFU
			if (fabs(ms.kx(i)/Ret) >= 2e-3 ||
				fabs(ms.kz(k)/Ret) >= 2e-2 ||
				(k==0 && i==0))
			{
				vo(2*i,jo1,k) = 0; vo(2*i+1,jo1,k) = 0;
				wo(2*i,jo1,k) = 0; wo(2*i+1,jo1,k) = 0;
				vo(2*i,jo2,k) = 0; vo(2*i+1,jo2,k) = 0;
				wo(2*i,jo2,k) = 0; wo(2*i+1,jo2,k) = 0;
			}
		}}

		uo.ifftx (jo1, jo2);
		vo.ifftxz(jo1, jo2);
		wo.ifftxz(jo1, jo2);
	}

	Bcond::SetBoundaryX(velo);
	Bcond::SetBoundaryZ(velo);

	delete[] hr;
	delete[] hi;

	// filter from MFU matching inner scale
	PIO::ReadSynthe(velb, time, Ret);

#ifdef MODULATION
	// modulation
	for (int j=0; j<=ms.Ny; j+=(j ? ms.Ny-1 : 1)) {

		double um = ub.meanxz(j);
		double vm = vb.meanxz(j);
		double wm = wb.meanxz(j);
		double uom= uo.meanxz(j<=1 ? jo1 : jo2);

		#pragma omp parallel for collapse(2)
		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			double x = ms.x(i), xc = ms.xc(i);
			double z = ms.z(k), zc = ms.zc(k);
			double y = j<=1 ? ms.yc(jo1) : ms.yc(jo2);

			double modu;
			double modv = gmav * (Interp::InterpNodeU(xc+dxv,y,zc+dzv,uo) - uom);
			double modw = gmaw * (Interp::InterpNodeU(xc+dxw,y,z +dzw,uo) - uom);

			if (ub(i,j,k) > um) modu = gmaup * (Interp::InterpNodeU(x+dxup,y,zc+dzup,uo) - uom);
			else                modu = gmaum * (Interp::InterpNodeU(x+dxum,y,zc+dzum,uo) - uom);

			if (j != 1) ub(i,j,k) = um + (ub(i,j,k) - um) * (1 + modu);
			if (j != 0) vb(i,j,k) = vm + (vb(i,j,k) - vm) * (1 + modv);
			if (j != 1) wb(i,j,k) = wm + (wb(i,j,k) - wm) * (1 + modw);
		}}
	}
#endif

	// combine small- & large-scales
	#pragma omp parallel for collapse(2)
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		ub(i,0,    k) += uo(i,jo1,k);
		ub(i,ms.Ny,k) += uo(i,jo2,k);

		double x = ms.xc(i) + dxav;
		double z = ms.zc(k) + dzav;

		vb(i,1,    k) += alfv * Interp::InterpNodeV(x,ms.y(jo1),z,vo);
		vb(i,ms.Ny,k) += alfv * Interp::InterpNodeV(x,ms.y(jo2),z,vo);

		x = ms.xc(i) + dxaw;
		z = ms.z (k) + dzaw;

		wb(i,0,    k) += alfw * Interp::InterpNodeW(x,ms.yc(jo1),z,wo);
		wb(i,ms.Ny,k) += alfw * Interp::InterpNodeW(x,ms.yc(jo2),z,wo);
	}}

	Bcond::SetBoundaryX(velb);
	Bcond::SetBoundaryZ(velb);
}


// Vctr PIO::Predict(double y, const Vctr &velout, const Vctr &veluni, double Ret, double rsclx, double rsclu, int cnt)
// {
// 	const Mesh &ms = velout.ms;
// 	const Mesh &ms0 = veluni.ms;

// 	double yb1 = y,               yb2 = 2. - y; // positions to predict
// 	double yo1 = 3.9 / sqrt(Ret), yo2 = 2. - yo1; // positions of outer signal

// 	int nx = (int)round(rsclx * ms.Lx / ms0.Lx) * (ms0.Nx-1) + 1;
// 	int nz = (int)round(rsclx * ms.Lz / ms0.Lz) * (ms0.Nz-1) + 1;

// 	Geometry_prdxz geo1(nx,   2,   nz, ms.Lx, yb2-yb1, ms.Lz); // to hold uS
// 	Geometry_prdxz geo2(ms.Nx,2,ms.Nz, ms.Lx, yo2-yo1, ms.Lz); // to hold uL

// 	geo1.InitMesh(0); geo1.InitInterval(); geo1.InitWaveNumber(); geo1.InitIndices();
// 	geo2.InitMesh(0); geo2.InitInterval(); geo2.InitWaveNumber(); geo2.InitIndices();

// 	const Mesh mss(geo1); Vctr vels(mss); vels.Set(0);
// 	const Mesh mso(geo2); Vctr velo(mso); velo.Set(0);

// 	Scla &us = vels[1], &vs = vels[2], &ws = vels[3];
// 	Scla &uo = velo[1], &vo = velo[2], &wo = velo[3];

// 	// PIO coefficients
// 	double *hr = new double[mso.Nxc], alfv, dxav, dzav;
// 	double *hi = new double[mso.Nxc], alfw, dxaw, dzaw;

// 	double gmaup, dxup, dzup;
// 	double gmaum, dxum, dzum;
// 	double gmav,  dxv,  dzv;
// 	double gmaw,  dxw,  dzw;

// 	read_PIO_Sup(hr, hi, alfv, dxav, dzav,
// 					 alfw, dxaw, dzaw, 1./Ret, yb1, geo2.kx, mso.Nxc);

// 	read_PIO_Mod(gmaup, dxup, dzup,
// 				 gmaum, dxum, dzum,
// 				 gmav,  dxv,  dzv,
// 				 gmaw,  dxw,  dzw, 1./Ret, yb1);

// 	// interpolate from MFU matching inner scale
// 	#pragma omp parallel for collapse(2)
// 	for (int j=0; j<=mss.Ny; j++) {
// 	for (int k=0; k<=mss.Nz; k++) {
// 	for (int i=0; i<=mss.Nx; i++) {

// 		double x = mss.x(i) * rsclx, xc = mss.xc(i) * rsclx;
// 		double z = mss.z(k) * rsclx, zc = mss.zc(k) * rsclx;
// 		double y = Filter::WallRscl(mss.y (j), rsclx);
// 		double yc= Filter::WallRscl(mss.yc(j), rsclx);

// 		if (i>0) us(i,j,k) = Interp::InterpNodeU(x,yc,zc,veluni[1]) * rsclu;
// 		if (j>0) vs(i,j,k) = Interp::InterpNodeV(xc,y,zc,veluni[2]) * rsclu;
// 		if (k>0) ws(i,j,k) = Interp::InterpNodeW(xc,yc,z,veluni[3]) * rsclu;

// 		// handle HALFMFU
// 		if (1+mss.Ny < 2/ms0.Ly*j) vs(i,j,k) *= -1;
// 	}}}

// 	// interpolate from LES outer region
// 	Interp::InterpBulkU(uo, velout[1]); // mean values will be removed later
// 	Interp::InterpBulkV(vo, velout[2]);
// 	Interp::InterpBulkW(wo, velout[3]);

// 	// get uL, vOL, wOL
// 	uo.fftx();
// 	vo.fftxz();
// 	wo.fftxz();

// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 		for (int j=0; j<=2; j+=2) {
// 		for (int i=0; i<mso.Nxc; i++) {

// 			// get mean value at wave number k_x = 0
// 			double um = 0;
// 			if (i == 0)
// 				for (int k=1; k<mso.Nz; k++)
// 					um += uo(2*i,j,k) / (mso.Nz-1.);

// 			for (int k=1; k<mso.Nz; k++) {
// 				double ur = uo(2*i,  j,k) - um;
// 				double ui = uo(2*i+1,j,k);
// 				// H needs conj because the fft defined in calibration is inversed
// 				uo(2*i,  j,k) = hr[i] * ur + hi[i] * ui;
// 				uo(2*i+1,j,k) = hr[i] * ui - hi[i] * ur;
// 			}
// 		}}

// 		#pragma omp for
// 		for (int k=0; k<mso.Nz-1; k++) {
// 		for (int i=0; i<mso.Nxc; i++) {
// 			// keep only the fluctuations with scale strictly larger than MFU
// 			if (fabs(mso.kx(i)/rsclx) >= 2.*PI/ms0.Lx ||
// 				fabs(mso.kz(k)/rsclx) >= 2.*PI/ms0.Lz ||
// 				(k==0 && i==0))
// 			{
// 				vo(2*i,1,k) = 0; vo(2*i+1,1,k) = 0;
// 				wo(2*i,0,k) = 0; wo(2*i+1,0,k) = 0;
// 				vo(2*i,2,k) = 0; vo(2*i+1,2,k) = 0;
// 				wo(2*i,2,k) = 0; wo(2*i+1,2,k) = 0;
// 			}
// 		}}
// 	}

// 	uo.ifftx();
// 	vo.ifftxz();
// 	wo.ifftxz();

// 	Bcond::SetBoundaryX(velo);
// 	Bcond::SetBoundaryZ(velo);

// 	delete[] hr;
// 	delete[] hi;

// 	// here MODULATION should always stay since the purpose of this function is to obtain highly resolved predicted field
// 	for (int j=0; j<=2; j++) {

// 		double um = us.meanxz(j);
// 		double vm = vs.meanxz(j);
// 		double wm = ws.meanxz(j);
// 		double uom = uo.meanxz(j); // uom seems unnecessary since mean of uo has been removed when calculating uL

// 		#pragma omp parallel for
// 		for (int k=1; k<mss.Nz; k++) {
// 		for (int i=1; i<mss.Nx; i++) {

// 			double x = mss.x(i), xc = mss.xc(i);
// 			double z = mss.z(k), zc = mss.zc(k);
// 			double y = j<=1 ? yo1 : yo2;

// 			double modu;
// 			double modv = gmav * (Interp::InterpNodeU(xc+dxv,y,zc+dzv,uo) - uom);
// 			double modw = gmaw * (Interp::InterpNodeU(xc+dxw,y,z +dzw,uo) - uom);

// 			if (us(i,j,k) > um) modu = gmaup * (Interp::InterpNodeU(x+dxup,y,zc+dzup,uo) - uom);
// 			else                modu = gmaum * (Interp::InterpNodeU(x+dxum,y,zc+dzum,uo) - uom);

// 			if (j != 1) us(i,j,k) = um + (us(i,j,k) - um) * (1 + modu);
// 			if (j != 0) vs(i,j,k) = vm + (vs(i,j,k) - vm) * (1 + modv);
// 			if (j != 1) ws(i,j,k) = wm + (ws(i,j,k) - wm) * (1 + modw);
// 		}}
// 	}

// 	// superposition
// 	#pragma omp parallel for
// 	for (int k=1; k<mss.Nz; k++) {
// 	for (int i=1; i<mss.Nx; i++) {

// 		double x = mss.x (i);
// 		double z = mss.zc(k);

// 		us(i,0,k) += Interp::InterpNodeU(x,yo1,z,uo);
// 		us(i,2,k) += Interp::InterpNodeU(x,yo2,z,uo);

// 		x = mss.xc(i) + dxav;
// 		z = mss.zc(k) + dzav;

// 		vs(i,1,k) += alfv * Interp::InterpNodeV(x,yo1,z,vo);
// 		vs(i,2,k) += alfv * Interp::InterpNodeV(x,yo2,z,vo);

// 		x = mss.xc(i) + dxaw;
// 		z = mss.z (k) + dzaw;

// 		ws(i,0,k) += alfw * Interp::InterpNodeW(x,yo1,z,wo);
// 		ws(i,2,k) += alfw * Interp::InterpNodeW(x,yo2,z,wo);
// 	}}

// 	Bcond::SetBoundaryX(vels);
// 	Bcond::SetBoundaryZ(vels);

// 	// record for debug
// 	char str[32];
// 	sprintf(str, "UBOT%08i", max(cnt,0)); us.debug_AsciiOutput("test/probedata/", str, 0,1);
// 	sprintf(str, "VBOT%08i", max(cnt,0)); vs.debug_AsciiOutput("test/probedata/", str, 1,2);
// 	sprintf(str, "WBOT%08i", max(cnt,0)); ws.debug_AsciiOutput("test/probedata/", str, 0,1);
// 	sprintf(str, "UTOP%08i", max(cnt,0)); us.debug_AsciiOutput("test/probedata/", str, 2,3);
// 	sprintf(str, "VTOP%08i", max(cnt,0)); vs.debug_AsciiOutput("test/probedata/", str, 2,3);
// 	sprintf(str, "WTOP%08i", max(cnt,0)); ws.debug_AsciiOutput("test/probedata/", str, 2,3);

// 	return vels;
// }


void PIO::PredictBoundarySGS(Vctr &shearsgs, const Vctr &velout, double Ret)
// add modulation effect to fluctuating SGS stress (since SGS is in itself small-scale, there is no superposition)
{
#ifdef MODULATION
	const Mesh &ms = shearsgs.ms;

	double yo1 = 3.9 / sqrt(Ret); // position of outer signal
	double yo2 = 2. - yo1;

	Geometry_prdxz geo(ms.Nx,2,ms.Nz, ms.Lx, yo2-yo1, ms.Lz); // to hold uL
	const Mesh mso(geo.Init(0));
	Vctr velo(mso);
	velo.Set(0);

	Scla &uo = velo[1];

	Scla &tau12sgs = shearsgs[1];
	Scla &tau23sgs = shearsgs[2];
	// Scla &tau13sgs = shearsgs[3];

	// interpolate from LES outer region
	Interp::InterpBulkU(uo, velout[1]); // mean values will be removed later

	// PIO coefficients
	double *hr = new double[mso.Nxc], alfv, dxav, dzav;
	double *hi = new double[mso.Nxc], alfw, dxaw, dzaw;

	double gmaup, dxup, dzup;
	double gmaum, dxum, dzum;
	double gmav,  dxv,  dzv;
	double gmaw,  dxw,  dzw;

	read_PIO_Sup(hr, hi, alfv, dxav, dzav,
					 alfw, dxaw, dzaw, 1./Ret, ms.y(1), geo.kx, mso.Nxc);

	read_PIO_Mod(gmaup, dxup, dzup,
				 gmaum, dxum, dzum,
				 gmav,  dxv,  dzv,
				 gmaw,  dxw,  dzw, 1./Ret, ms.y(1));

	// get uL, vOL, wOL
	uo.fftx();

	#pragma omp parallel for collapse(2)
	for (int j=0; j<=2; j+=2) {
	for (int i=0; i<mso.Nxc; i++) {
		// get mean value at wave number k_x = 0
		double um = 0;
		if (i == 0)
			for (int k=1; k<mso.Nz; k++)
				um += uo(2*i,j,k) / (mso.Nz-1.);

		for (int k=1; k<mso.Nz; k++) {
			double ur = uo(2*i,  j,k) - um;
			double ui = uo(2*i+1,j,k);
			// H needs conj because the fft defined in calibration is inversed
			uo(2*i,  j,k) = hr[i] * ur + hi[i] * ui;
			uo(2*i+1,j,k) = hr[i] * ui - hi[i] * ur;
		}
	}}

	uo.ifftx();

	Bcond::SetBoundaryX(velo);
	Bcond::SetBoundaryZ(velo);

	// modulation
	#pragma omp parallel for collapse(2)
	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		double x = ms.x(i), xc = ms.xc(i);
		double z = ms.z(k), zc = ms.zc(k);
		double y = j<=1 ? yo1 : yo2;

		double modw   = gmaw * Interp::InterpNodeU(xc+dxw,y,z+dzw,uo);
		double modv23 = gmav * Interp::InterpNodeU(xc+dxv,y,z+dzv,uo);
		double modv12 = gmav * Interp::InterpNodeU(x+dxv,y,zc+dzv,uo);
		double modup = gmaup * Interp::InterpNodeU(x+dxup,y,zc+dzup,uo);
		double modum = gmaum * Interp::InterpNodeU(x+dxum,y,zc+dzum,uo);
		double modu = .5 * (modup + modum);

		if (i > 0) tau12sgs(i,j,k) *= 1 + modu + modv12;
		if (k > 0) tau23sgs(i,j,k) *= 1 + modw + modv23;
	}}}

	delete[] hr;
	delete[] hi;
#endif
}



void read_PIO_Sup(double *hr, double *hi,
	double &alfv, double &dxav, double &dzav,
	double &alfw, double &dxaw, double &dzaw,
	double lc, double y, const double *kx, int Nxc)
// get PIO coefficients
{
	// read files
	char str[1024];

	int nx, ny;

	FILE *fp = fopen("calibration_PIO/U/HL.dat", "r");
	fgets(str, 1024, fp);
	fgets(str, 1024, fp);
	fgets(str, 1024, fp);

	sscanf(strstr(str, "i ="), "i = %i", &nx);
	sscanf(strstr(str, "j ="), "j = %i", &ny);

	double *kxp = new double[nx];
	double *ysp = new double[ny];
	double *hrs = new double[nx*ny];
	double *his = new double[nx*ny];

	for (int j=0; j<ny; j++) {
	for (int i=0; i<nx; i++) {
		fgets(str, 1024, fp);
		sscanf(str, "%le %le %le %le",
			&kxp[i], &ysp[j], &hrs[j*nx+i], &his[j*nx+i]);
	}}
	fclose(fp);

	double *alfvs = new double[ny];
	double *dxbvs = new double[ny];
	double *dzbvs = new double[ny];
	double *alfws = new double[ny];
	double *dxbws = new double[ny];
	double *dzbws = new double[ny];

	fp = fopen("calibration_PIO/V/Alpha.dat", "r");
	for (int j=-3; j<ny; j++) {
		fgets(str, 1024, fp);
		if (j >= 0)
			sscanf(str, "%le %le %le %le",
				&ysp[j], &alfvs[j], &dxbvs[j], &dzbvs[j]);
	}
	fclose(fp);

	fp = fopen("calibration_PIO/W/Alpha.dat", "r");
	for (int j=-3; j<ny; j++) {
		fgets(str, 1024, fp);
		if (j >= 0)
			sscanf(str, "%le %le %le %le",
				&ysp[j], &alfws[j], &dxbws[j], &dzbws[j]);
	}
	fclose(fp);

	// interpolation
	double yp = y / lc;
	int j0 = Interp::BiSearch(yp, ysp, 0, ny-1);

	double a = ysp[j0+1] - yp;
	double b = yp - ysp[j0];

	alfv = (a * alfvs[j0] + b * alfvs[j0+1]) / (a+b);
	dxav = (a * dxbvs[j0] + b * dxbvs[j0+1]) / (a+b) * lc;
	dzav = (a * dzbvs[j0] + b * dzbvs[j0+1]) / (a+b) * lc;

	alfw = (a * alfws[j0] + b * alfws[j0+1]) / (a+b);
	dxaw = (a * dxbws[j0] + b * dxbws[j0+1]) / (a+b) * lc;
	dzaw = (a * dzbws[j0] + b * dzbws[j0+1]) / (a+b) * lc;


	for (int i=0; i<Nxc; i++) {

		double kp = fabs(kx[i]) * lc;
		int i0 = Interp::BiSearch(kp, kxp, 0, nx-1);

		double c = kxp[i0+1] - kp;
		double d = kp - kxp[i0];

		hr[i] = 1. / ((a+b) * (c+d)) * (
			a*c * hrs[ j0   *nx + i0  ] +
			a*d * hrs[ j0   *nx + i0+1] +
			b*c * hrs[(j0+1)*nx + i0  ] +
			b*d * hrs[(j0+1)*nx + i0+1] );

		hi[i] = 1. / ((a+b) * (c+d)) * (
			a*c * his[ j0   *nx + i0  ] +
			a*d * his[ j0   *nx + i0+1] +
			b*c * his[(j0+1)*nx + i0  ] +
			b*d * his[(j0+1)*nx + i0+1] );
	}

	delete[] kxp;
	delete[] ysp;
	delete[] hrs;
	delete[] his;
	delete[] alfvs; delete[] dxbvs; delete[] dzbvs;
	delete[] alfws; delete[] dxbws; delete[] dzbws;
}



void read_PIO_Mod(
	double &gmaup, double &dxup, double &dzup,
	double &gmaum, double &dxum, double &dzum,
	double &gmav,  double &dxv,  double &dzv,
	double &gmaw,  double &dxw,  double &dzw,
	double lc, double y)
// get PIO modulation coefficients
{
	// read files
	char str[1024];

	int ny;

	FILE *fp = fopen("calibration_PIO/U/UP/Gamma.dat", "r");
	fgets(str, 1024, fp);
	fgets(str, 1024, fp);
	fgets(str, 1024, fp);

	sscanf(strstr(str, "i ="), "i = %i", &ny);

	double *ysp = new double[ny];
	double *gmaups = new double[ny], *dxups = new double[ny], *dzups = new double[ny];
	double *gmaums = new double[ny], *dxums = new double[ny], *dzums = new double[ny];
	double *gmavs  = new double[ny], *dxvs  = new double[ny], *dzvs  = new double[ny];
	double *gmaws  = new double[ny], *dxws  = new double[ny], *dzws  = new double[ny];

	for (int j=0; j<ny; j++) {
		fgets(str, 1024, fp);
		sscanf(str, "%le %le %le %le",
			&ysp[j], &gmaups[j], &dxups[j], &dzups[j]);
	}
	fclose(fp);

	fp = fopen("calibration_PIO/U/UM/Gamma.dat", "r");
	for (int j=-3; j<ny; j++) {
		fgets(str, 1024, fp);
		if (j >= 0)
			sscanf(str, "%le %le %le %le",
				&ysp[j], &gmaums[j], &dxums[j], &dzums[j]);
	}
	fclose(fp);


	fp = fopen("calibration_PIO/V/Gamma.dat", "r");
	for (int j=-3; j<ny; j++) {
		fgets(str, 1024, fp);
		if (j >= 0)
			sscanf(str, "%le %le %le %le",
				&ysp[j], &gmavs[j], &dxvs[j], &dzvs[j]);
	}
	fclose(fp);

	fp = fopen("calibration_PIO/W/Gamma.dat", "r");
	for (int j=-3; j<ny; j++) {
		fgets(str, 1024, fp);
		if (j >= 0)
			sscanf(str, "%le %le %le %le",
				&ysp[j], &gmaws[j], &dxws[j], &dzws[j]);
	}
	fclose(fp);

	// interpolation
	double yp = y / lc;
	int j0 = Interp::BiSearch(yp, ysp, 0, ny-1);

	double a = ysp[j0+1] - yp;
	double b = yp - ysp[j0];

	gmaup = (a * gmaups[j0] + b * gmaups[j0+1]) / (a+b);
	dxup  = (a * dxups [j0] + b * dxups [j0+1]) / (a+b) * lc;
	dzup  = (a * dzups [j0] + b * dzups [j0+1]) / (a+b) * lc;

	gmaum = (a * gmaums[j0] + b * gmaums[j0+1]) / (a+b);
	dxum  = (a * dxums [j0] + b * dxums [j0+1]) / (a+b) * lc;
	dzum  = (a * dzums [j0] + b * dzums [j0+1]) / (a+b) * lc;

	gmav = (a * gmavs[j0] + b * gmavs[j0+1]) / (a+b);
	dxv  = (a * dxvs [j0] + b * dxvs [j0+1]) / (a+b) * lc;
	dzv  = (a * dzvs [j0] + b * dzvs [j0+1]) / (a+b) * lc;

	gmaw = (a * gmaws[j0] + b * gmaws[j0+1]) / (a+b);
	dxw  = (a * dxws [j0] + b * dxws [j0+1]) / (a+b) * lc;
	dzw  = (a * dzws [j0] + b * dzws [j0+1]) / (a+b) * lc;


	delete[] gmaups; delete[] dxups; delete[] dzups;
	delete[] gmaums; delete[] dxums; delete[] dzums;
	delete[] gmavs;  delete[] dxvs;  delete[] dzvs;
	delete[] gmaws;  delete[] dxws;  delete[] dzws;
}


void PIO::ReadSynthe(Vctr &vel, double time, double Ret)
// read synthetic field from '.bin' files to boundaries of vel
{
	// read info section
	int   n1, n2, n3;
	float dx, dz, dt;

	ifstream us("probedata/us.bin", ios::in|ios::binary);
	ifstream vs("probedata/vs.bin", ios::in|ios::binary);
	ifstream ws("probedata/ws.bin", ios::in|ios::binary);

	us.read((char*)(&n1), sizeof(int));
	us.read((char*)(&n2), sizeof(int));
	us.read((char*)(&n3), sizeof(int));
	us.read((char*)(&dx), sizeof(float));
	us.read((char*)(&dz), sizeof(float));
	us.read((char*)(&dt), sizeof(float));

	dx /= Ret;
	dz /= Ret;
	dt /= Ret;

	// prepare space to hold data
	const Mesh &ms = vel.ms;
	Scla &u = vel[1];
	Scla &v = vel[2];
	Scla &w = vel[3];

	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		u(i,0,k) = u(i,ms.Ny,k) = 0;
		v(i,1,k) = v(i,ms.Ny,k) = 0;
		w(i,0,k) = w(i,ms.Ny,k) = 0;
	}}

	// file data input and interpolate to velocity
	double *buf = new double[n1*n2];
	double *bvf = new double[n1*n2];
	double *bwf = new double[n1*n2];

	int    t1 = fmod(time/dt,n3), t2 = (t1 + 1) % n3;
	double c1 = fmod(time,dt)/dt, c2 = 1. - c1;

	// double r12 = 0;
	
	for (int t=t1; t>=0; t=(t==t2 ? -1 : t2)) {

		double c = t==t2 ? c1 : c2;

		us.seekg(n1*n2*(t+1) * sizeof(double), ios::beg);
		vs.seekg(n1*n2*(t+1) * sizeof(double), ios::beg);
		ws.seekg(n1*n2*(t+1) * sizeof(double), ios::beg);

		us.read((char*)buf, n1*n2 * sizeof(double));
		vs.read((char*)bvf, n1*n2 * sizeof(double));
		ws.read((char*)bwf, n1*n2 * sizeof(double));

		// for (int k1=0; k1<n2; k1++) {
		// for (int i1=0; i1<n1; i1++) {
		// 	r12 += c/(n1*n2) * (buf[k1*n1+i1] * bvf[k1*n1+i1]);
		// }}

		#pragma omp parallel for collapse(2)
		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {
			// interpolation must be precise to ensure correct Reynolds stress
			int i1 = Interp::ShiftPrd(ms.x(i), 0, n1*dx) / dx, i2 = (i1+1)%n1;
			int k1 = Interp::ShiftPrd(ms.z(k), 0, n2*dz) / dz, k2 = (k1+1)%n2;

			double a1 = fmod(fmod(ms.x(i), dx) + dx, dx) / dx, a2 = 1. - a1;
			double b1 = fmod(fmod(ms.z(k), dz) + dz, dz) / dz, b2 = 1. - b1;

			double ub = c * (
				b2 * (a2 * buf[k1*n1+i1] + a1 * buf[k1*n1+i2]) +
				b1 * (a2 * buf[k2*n1+i1] + a1 * buf[k2*n1+i2]) );

			double vb = c * (
				b2 * (a2 * bvf[k1*n1+i1] + a1 * bvf[k1*n1+i2]) +
				b1 * (a2 * bvf[k2*n1+i1] + a1 * bvf[k2*n1+i2]) );

			double wb = c * (
				b2 * (a2 * bwf[k1*n1+i1] + a1 * bwf[k1*n1+i2]) +
				b1 * (a2 * bwf[k2*n1+i1] + a1 * bwf[k2*n1+i2]) );

			u(i,0,k) += ub; u(i,ms.Ny,k) += ub;
			v(i,1,k) += vb; v(i,ms.Ny,k) -= vb;
			w(i,0,k) += wb; w(i,ms.Ny,k) += wb;
		}}
	}

	// FILE* fp = fopen("r12ref.dat", "a");
	// fprintf(fp, "%.6e\n", r12);
	// fclose(fp);
}







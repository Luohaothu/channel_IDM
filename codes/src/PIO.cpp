#include "PIO.h"
#include "Interp.h"
#include "Filter.h"
#include "Bcond.h"


using namespace std;


void read_PIO(double *hr, double *hi,
	double &alfv, double &dxav, double &dzav,
	double &alfw, double &dxaw, double &dzaw,
	double lc, double y, const double *kx, int Nxc);


Vctr PIO::BoundaryPredict(const Vctr &vel, const Vctr &velmfu, double Ret, double rsclx, double rsclu)
{
	const Mesh &ms = vel.ms;
	const Mesh &ms0 = velmfu.ms;

	double yb1 = ms.y(1); // position to predict
	double yb2 = ms.y(ms.Ny);

	double yo1 = 3.9 / sqrt(Ret); // position of outer signal
	double yo2 = 2. - yo1;

	int nx = (int)round(rsclx * ms.Lx / ms0.Lx) * (ms0.Nx-1) + 1;
	int nz = (int)round(rsclx * ms.Lz / ms0.Lz) * (ms0.Nz-1) + 1;

	Geometry_prdxz geo1(nx,   2,   nz, ms.Lx, yb2-yb1, ms.Lz); // to hold uS
	Geometry_prdxz geo2(ms.Nx,2,ms.Nz, ms.Lx, yo2-yo1, ms.Lz); // to hold uL
	Geometry_prdxz geo3(ms.Nx,2,ms.Nz, ms.Lx, yb2-yb1, ms.Lz); // to hold u_p

	geo1.InitMesh(0); geo1.InitInterval(); geo1.InitWaveNumber(); geo1.InitIndices();
	geo2.InitMesh(0); geo2.InitInterval(); geo2.InitWaveNumber(); geo2.InitIndices();
	geo3.InitMesh(0); geo3.InitInterval(); geo3.InitWaveNumber(); geo3.InitIndices();

	const Mesh mss(geo1); Vctr vels(mss); vels.Set(0);
	const Mesh mso(geo2); Vctr velo(mso); velo.Set(0);
	const Mesh msb(geo3); Vctr velb(msb); velb.Set(0);

	Scla &us = vels[1], &vs = vels[2], &ws = vels[3];
	Scla &uo = velo[1], &vo = velo[2], &wo = velo[3];
	Scla &ub = velb[1], &vb = velb[2], &wb = velb[3];

	// PIO coefficients
	double *hr = new double[mso.Nxc], alfv, dxav, dzav;
	double *hi = new double[mso.Nxc], alfw, dxaw, dzaw;

	read_PIO(hr, hi, alfv, dxav, dzav,
					 alfw, dxaw, dzaw, 1./Ret, yb1, geo2.kx, mso.Nxc);

	// interpolate from MFU matching inner scale
	#pragma omp parallel for collapse(2)
	for (int j=0; j<=mss.Ny; j++) {
	for (int k=0; k<=mss.Nz; k++) {
	for (int i=0; i<=mss.Nx; i++) {

		double x = mss.x(i) * rsclx, xc = mss.xc(i) * rsclx;
		double z = mss.z(k) * rsclx, zc = mss.zc(k) * rsclx;
		double y = WallRscl(mss.y (j), rsclx);
		double yc= WallRscl(mss.yc(j), rsclx);

		if (i>0) us(i,j,k) = Interp::InterpNodeU(x,yc,zc,velmfu[1]) * rsclu;
		if (j>0) vs(i,j,k) = Interp::InterpNodeV(xc,y,zc,velmfu[2]) * rsclu;
		if (k>0) ws(i,j,k) = Interp::InterpNodeW(xc,yc,z,velmfu[3]) * rsclu;
	}}}

	// interpolate from LES outer region
	Interp::InterpBulkU(uo, vel[1]); // mean values will be removed later
	Interp::InterpBulkV(vo, vel[2]);
	Interp::InterpBulkW(wo, vel[3]);

	// get uOL, vOL, wOL
	uo.fftx();
	vo.fftxz();
	wo.fftxz();

	#pragma omp parallel
	{
		#pragma omp for collapse(2)
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

		#pragma omp for
		for (int k=0; k<mso.Nz-1; k++) {
		for (int i=0; i<mso.Nxc; i++) {
			// keep only the fluctuations with scale strictly larger than MFU
			if (fabs(mso.kx(i)/rsclx) >= 2.*PI/ms0.Lx ||
				fabs(mso.kz(k)/rsclx) >= 2.*PI/ms0.Lz ||
				(k==0 && i==0))
			{
				vo(2*i,1,k) = 0; vo(2*i+1,1,k) = 0;
				wo(2*i,0,k) = 0; wo(2*i+1,0,k) = 0;
				vo(2*i,2,k) = 0; vo(2*i+1,2,k) = 0;
				wo(2*i,2,k) = 0; wo(2*i+1,2,k) = 0;
			}
		}}
	}

	uo.ifftx();
	vo.ifftxz();
	wo.ifftxz();

	Bcond::SetBoundaryX(velo);
	Bcond::SetBoundaryZ(velo);

	// combine small- & large-scales
	#pragma omp parallel for
	for (int k=1; k<msb.Nz; k++) {
	for (int i=1; i<msb.Nx; i++) {

		// small-scale component
		double x = msb.x(i), xc = msb.xc(i), dx = msb.dx(i), hx = msb.hx(i);
		double z = msb.z(k), zc = msb.zc(k), dz = msb.dz(k), hz = msb.hz(k);

		ub(i,0,k) = Filter::FilterNodeU(x, yb1,zc, hx,0,dz, us);
		vb(i,1,k) = Filter::FilterNodeV(xc,yb1,zc, dx,0,dz, vs);
		wb(i,0,k) = Filter::FilterNodeW(xc,yb1,z,  dx,0,hz, ws);
		ub(i,2,k) = Filter::FilterNodeU(x, yb2,zc, hx,0,dz, us);
		vb(i,2,k) = Filter::FilterNodeV(xc,yb2,zc, dx,0,dz, vs);
		wb(i,2,k) = Filter::FilterNodeW(xc,yb2,z,  dx,0,hz, ws);

		// large-scale component
		ub(i,0,k) += uo(i,0,k);
		ub(i,2,k) += uo(i,2,k);

		x = msb.xc(i) + dxav;
		z = msb.zc(k) + dzav;

		vb(i,1,k) += alfv * Interp::InterpNodeV(x,yo1,z,vo);
		vb(i,2,k) += alfv * Interp::InterpNodeV(x,yo2,z,vo);

		x = msb.xc(i) + dxaw;
		z = msb.z (k) + dzaw;

		wb(i,0,k) += alfw * Interp::InterpNodeW(x,yo1,z,wo);
		wb(i,2,k) += alfw * Interp::InterpNodeW(x,yo2,z,wo);
	}}

	Bcond::SetBoundaryX(velb);
	Bcond::SetBoundaryZ(velb);

	delete[] hr;
	delete[] hi;

	// record for debug
	static int cnt = 0;

	if ((++cnt) % 1000 == 0) {
	
		#pragma omp parallel for
		for (int k=1; k<mss.Nz; k++) {
		for (int i=1; i<mss.Nx; i++) {

			double x = mss.x (i);
			double z = mss.zc(k);

			us(i,0,k) += Interp::InterpNodeU(x,yo1,z,uo);
			us(i,2,k) += Interp::InterpNodeU(x,yo2,z,uo);

			x = mss.xc(i) + dxav;
			z = mss.zc(k) + dzav;

			vs(i,1,k) += alfv * Interp::InterpNodeV(x,yo1,z,vo);
			vs(i,2,k) += alfv * Interp::InterpNodeV(x,yo2,z,vo);

			x = mss.xc(i) + dxaw;
			z = mss.z (k) + dzaw;

			ws(i,0,k) += alfw * Interp::InterpNodeW(x,yo1,z,wo);
			ws(i,2,k) += alfw * Interp::InterpNodeW(x,yo2,z,wo);
		}}

		Bcond::SetBoundaryX(vels);
		Bcond::SetBoundaryZ(vels);

		char str[32];

		sprintf(str, "UBOT%08i", cnt); us.debug_AsciiOutput("test/probedata/", str, 0,1);
		sprintf(str, "VBOT%08i", cnt); vs.debug_AsciiOutput("test/probedata/", str, 1,2);
		sprintf(str, "WBOT%08i", cnt); ws.debug_AsciiOutput("test/probedata/", str, 0,1);
		sprintf(str, "UTOP%08i", cnt); us.debug_AsciiOutput("test/probedata/", str, 2,3);
		sprintf(str, "VTOP%08i", cnt); vs.debug_AsciiOutput("test/probedata/", str, 2,3);
		sprintf(str, "WTOP%08i", cnt); ws.debug_AsciiOutput("test/probedata/", str, 2,3);
	}

	return velb;
}





void read_PIO(double *hr, double *hi,
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
	double yp = WallRscl(y, 1./lc);
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








# include <iostream>
# include <stdlib.h>
# include <cmath>

# include "Statis.h"

using namespace std;


Statis::Statis(const Mesh &ms):
ms(ms)
{	
	um_ = new double[ms.Ny+1];
	vm_ = new double[ms.Ny+1];
	wm_ = new double[ms.Ny+1];
	pm_ = new double[ms.Ny+1];
	r11_= new double[ms.Ny+1];
	r22_= new double[ms.Ny+1];
	r33_= new double[ms.Ny+1];
	r12_= new double[ms.Ny+1];
	r23_= new double[ms.Ny+1];
	r13_= new double[ms.Ny+1];
	rpu_= new double[ms.Ny+1];
	rpv_= new double[ms.Ny+1];
	rpw_= new double[ms.Ny+1];
	rpp_= new double[ms.Ny+1];
	num_= new double[ms.Ny+1];
}

Statis::~Statis()
{
	delete[] um_;  delete[] vm_;  delete[] wm_;
	delete[] r11_; delete[] r22_; delete[] r33_;
	delete[] r12_; delete[] r23_; delete[] r13_;
	delete[] rpu_; delete[] rpv_; delete[] rpw_;
	delete[] pm_;  delete[] rpp_; delete[] num_;
}


double Statis::ReynoldsStress(int j, const Vctr &vel)
/* compute the Reynolds stress R12 on the edge for a channel flow */
{
	const Mesh &ms = vel.ms;
	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	double rs = 0;
	double um = 0, ua;
	double vm = 0, va;

	#pragma omp parallel for reduction(+: rs,um,vm) collapse(2)
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		ua = .5/ms.hy(j) * (u(i,j,k) * ms.dy(j-1) + u(i,ms.jma(j),k) * ms.dy(j));
		va = .5/ms.hx(i) * (v(i,j,k) * ms.dx(i-1) + v(ms.ima(i),j,k) * ms.dx(i));

		rs += ua * va;
		um += ua;
		vm += va;
	}}

	double weight = 1. / (ms.Nx-1.) / (ms.Nz-1.);

	rs -= weight * um * vm;
	rs *= weight;

	return rs;
}

vector<double> Statis::ReynoldsShearStresses(int j, const Vctr &vel)
/* compute three Reynolds shear stresses on the edges for a channel flow */
{
	const Mesh &ms = vel.ms;
	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	double rs12=0, um12=0, vm12=0, wa;
	double rs23=0, vm23=0, wm23=0, ua;
	double rs13=0, um13=0, wm13=0, va;

	#pragma omp parallel for reduction(+: rs12,rs23,rs13,um12,vm12,vm23,wm23,um13,wm13) collapse(2)
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		// interpolate to Z-edge
		ua = .5/ms.hy(j) * (u(i,j,k) * ms.dy(j-1) + u(i,ms.jma(j),k) * ms.dy(j));
		va = .5/ms.hx(i) * (v(i,j,k) * ms.dx(i-1) + v(ms.ima(i),j,k) * ms.dx(i));

		rs12 += ua * va;
		um12 += ua;
		vm12 += va;

		// interpolate to X-edge
		va = .5/ms.hz(k) * (v(i,j,k) * ms.dz(k-1) + v(i,j,ms.kma(k)) * ms.dz(k));
		wa = .5/ms.hy(j) * (w(i,j,k) * ms.dy(j-1) + w(i,ms.jma(j),k) * ms.dy(j));

		rs23 += va * wa;
		vm23 += ua;
		wm23 += va;

		// interpolate to Y_edge
		ua = .5/ms.hz(k) * (u(i,j,k) * ms.dz(k-1) + u(i,j,ms.kma(k)) * ms.dz(k));
		wa = .5/ms.hx(i) * (w(i,j,k) * ms.dx(i-1) + w(ms.ima(i),j,k) * ms.dx(i));

		rs13 += ua * wa;
		um13 += ua;
		wm13 += wa;
	}}

	double weight = 1. / (ms.Nx-1.) / (ms.Nz-1.);

	rs12 -= weight * um12 * vm12;  rs12 *= weight;
	rs23 -= weight * vm23 * wm23;  rs23 *= weight;
	rs13 -= weight * um13 * wm13;  rs13 *= weight;

	vector<double> rs = {rs12, rs23, rs13};

	return rs;
}

const double* Statis::ReynoldsNormalStresses(int j, const Vctr &vel)
/* compute three Reynolds normal stresses on cell centers for a channel flow */
{
	const Mesh &ms = vel.ms;
	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	double rs11=0, um=0, ua;
	double rs22=0, vm=0, va;
	double rs33=0, wm=0, wa;

	#pragma omp parallel for reduction(+: rs11,rs22,rs33,um,vm,wm) collapse(2)
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {
		// interpolate to cell centers
		ua = .5 * (u(i,j,k) + u(ms.ipa(i),j,k));
		va = .5 * (v(i,j,k) + v(i,ms.jpa(j),k));
		wa = .5 * (w(i,j,k) + w(i,j,ms.kpa(k)));

		rs11 += ua * ua;  um += ua;
		rs22 += va * va;  vm += va;
		rs33 += wa * wa;  wm += wa;
	}}

	double weight = 1. / (ms.Nx-1.) / (ms.Nz-1.);

	rs11 -= weight * um * um;  rs11 *= weight;
	rs22 -= weight * vm * vm;  rs22 *= weight;
	rs33 -= weight * wm * wm;  rs33 *= weight;

	static double rs[3];
	rs[0] = rs11;
	rs[1] = rs22;
	rs[2] = rs33;

	return rs;
}

const double* Statis::MeanVisShearStresses(int j, const Vctr &vel, const Flow &vis)
// compute mean viscous (including kinetic and eddy viscosity) shear stresses on edges for a channel flow
{
	const Mesh &ms = vel.ms;

	const Scla &u = vel[1], &nux = vis.SeeVec(1);
	const Scla &v = vel[2], &nuy = vis.SeeVec(2);
	const Scla &w = vel[3], &nuz = vis.SeeVec(3);

	double tau12 = 0;
	double tau23 = 0;
	double tau13 = 0;

	#pragma omp parallel for reduction(+: tau12,tau23,tau13) collapse(2)
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		const double *sr = vel.ShearStrain(i,j,k);

		tau12 += 2. * nuz(i,j,k) * sr[0];
		tau23 += 2. * nux(i,j,k) * sr[1];
		tau13 += 2. * nuy(i,j,k) * sr[2];
	}}

	double weight = 1. / (ms.Nx-1.) / (ms.Nz-1.);

	static double tau[3];

	tau[0] = tau12 * weight;
	tau[1] = tau23 * weight;
	tau[2] = tau13 * weight;

	return tau;
}

const double* Statis::MeanVisNormalStresses(int j, const Vctr &vel, const Flow &vis)
// compute mean viscous (including kinetic and eddy viscosity) normal stresses on cell centers for a channel flow
{
	const Mesh &ms = vel.ms;

	const Scla &u = vel[1];
	const Scla &v = vel[2];
	const Scla &w = vel[3];

	double tau11 = 0;
	double tau22 = 0;
	double tau33 = 0;

	#pragma omp parallel for reduction(+: tau11,tau22,tau33) collapse(2)
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		double nu = vis.SeeScl()(i,j,k);

		tau11 += 2. * nu * (u(ms.ipa(i),j,k) - u(i,j,k)) / ms.dx(i);
		tau22 += 2. * nu * (v(i,ms.jpa(j),k) - v(i,j,k)) / ms.dy(j);
		tau33 += 2. * nu * (w(i,j,ms.kpa(k)) - w(i,j,k)) / ms.dz(k);
	}}

	double weight = 1. / (ms.Nx-1.) / (ms.Nz-1.);

	static double tau[3];

	tau[0] = tau11 * weight;
	tau[1] = tau22 * weight;
	tau[2] = tau33 * weight;

	return tau;
}

void Statis::Check(const Flow &fld, const Flow &vis, double Re, double dt)
{
	div_ = CheckDiv (fld.SeeVec(), divpos_);
	cfl_ = CheckCfl (fld.SeeVec(), dt, cflpos_);
	CheckTaub(fld.SeeVec(), vis, taub_);
	ener_ = CheckProf(fld, vis, velm_);
}


double Statis::CheckDiv(const Vctr &vel, int pos[3])
{
	const Mesh &ms = vel.ms;

	double div, divmax = 0;

	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		if (divmax < (div = fabs(vel.Divergence(i,j,k)))) {

			divmax = div;
			pos[0] = i;
			pos[1] = j;
			pos[2] = k;
		}
	}}}
	return divmax;
}

double Statis::CheckCfl(const Vctr &vel, double dt, int pos[3])
{
	const Mesh &ms = vel.ms;

	double cfl, cflmax = 0;

	for (int j=1; j<ms.Ny; j++) {
	for (int k=1; k<ms.Nz; k++) {
	for (int i=1; i<ms.Nx; i++) {

		if (cflmax < (cfl = vel.Convection(i,j,k) * dt)) {

			cflmax = cfl;
			pos[0] = i;
			pos[1] = j;
			pos[2] = k;
		}
	}}}
	return cflmax;
}

void Statis::CheckTaub(const Vctr &vel, const Flow &vis, double taub[3])
/* calculate 3 components of total stress acting on real boundary:
   tau_2i = 2 * nu < du_i/dy + dv/dx_i > - <v'u_i'> */
{
	const Mesh &ms = vel.ms;

	double tau21 = 0;
	double tau22 = 0;
	double tau23 = 0;

	{ // lower boundary
		const vector<double> rss  = Statis::ReynoldsShearStresses (1, vel);
		const double *rsn  = Statis::ReynoldsNormalStresses(1, vel);
		const double *taus = Statis::MeanVisShearStresses  (1, vel, vis);
		const double *taun = Statis::MeanVisNormalStresses (1, vel, vis);

		tau21 += .5 * (taus[0] - rss[0]);
		tau22 += .5 * (taun[1] - rsn[1]); // use yc[1] to approximate stress at y[1]
		tau23 += .5 * (taus[1] - rss[1]);
	}{ // upper boundary
		const vector<double> rss  = Statis::ReynoldsShearStresses (ms.Ny,   vel);
		const double *rsn  = Statis::ReynoldsNormalStresses(ms.Ny-1, vel);
		const double *taus = Statis::MeanVisShearStresses  (ms.Ny,   vel, vis);
		const double *taun = Statis::MeanVisNormalStresses (ms.Ny-1, vel, vis);

		tau21 -= .5 * (taus[0] - rss[0]);
		tau22 += .5 * (taun[1] - rsn[1]); // use yc[Ny-1] to approximate stress at y[Ny]
		tau23 -= .5 * (taus[1] - rss[1]);
	}

	taub[0] = tau21;
	taub[1] = tau22;
	taub[2] = tau23;
}

double Statis::CheckProf(const Flow &fld, const Flow &vis, double velm[3])
/* calculate mean frofiles, Reynolds stresses and pressure-velocity correlations,
   aa well as bulk velocity and total fluctuation energy, at cell centers */
{
	const Mesh &ms = fld.ms;

	Scla u(ms);
	Scla v(ms);
	Scla w(ms);
	Scla p(ms);
	Scla nuc(ms);

	// interpolate to cell-centered points
	fld.SeeVec(1).Ugrid2CellCenter(u);
	fld.SeeVec(2).Vgrid2CellCenter(v);
	fld.SeeVec(3).Wgrid2CellCenter(w);
	
	p.Set(fld.SeeScl());
	p.SetLyr(p.SeeLyr(1), 0);
	p.SetLyr(p.SeeLyr(ms.Ny-1), ms.Ny);

	nuc.Set(vis.SeeScl());

	// reset accumulators
	double ener = 0;

	velm[0] = 0;
	velm[1] = 0;
	velm[2] = 0;

	for (int j=0; j<=ms.Ny; j++) {
		r11_[j]=0; r22_[j]=0; r33_[j]=0;
		r12_[j]=0; r23_[j]=0; r13_[j]=0;
		rpu_[j]=0; rpv_[j]=0; rpw_[j]=0;
		rpp_[j]=0;	
	}

	// calculate statistics at cell centers
	for (int j=0; j<=ms.Ny; j++) {

		// calculate means & fluctuations
		u.MnsLyr( um_[j] = u.meanxz(j), j );
		v.MnsLyr( vm_[j] = v.meanxz(j), j );
		w.MnsLyr( wm_[j] = w.meanxz(j), j );
		p.MnsLyr( pm_[j] = p.meanxz(j), j );

		num_[j] = nuc.meanxz(j);

		// calculate Reynolds stress & cross correlations
		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			int id = ms.idx(i,j,k);
			double weight = 1. / (ms.Nx-1) / (ms.Nz-1);
			
			r11_[j] += u[id] * u[id] * weight;
			r22_[j] += v[id] * v[id] * weight;
			r33_[j] += w[id] * w[id] * weight;
			r12_[j] += u[id] * v[id] * weight;
			r23_[j] += v[id] * w[id] * weight;
			r13_[j] += u[id] * w[id] * weight;
			rpu_[j] += p[id] * u[id] * weight;
			rpv_[j] += p[id] * v[id] * weight;
			rpw_[j] += p[id] * w[id] * weight;
			rpp_[j] += p[id] * p[id] * weight;
		}}

		// calculate bulk velocity & turbulent energy (wall-normal integration)
		double weight = (0<j && j<ms.Ny) ? ms.dy(j)/ms.Ly : 0;

		velm[0] += um_[j] * weight;
		velm[1] += vm_[j] * weight;
		velm[2] += wm_[j] * weight;

		ener += .5*weight * (r11_[j] + r22_[j] + r33_[j]);
	}

	return ener;
}


void Statis::WriteProfile(const char *path, int tstep) const
{
	char str[1024];

	if (tstep < 0)	sprintf(str, "%sPROF.dat",     path?path:"");
	else			sprintf(str, "%sPROF%08i.dat", path?path:"", tstep);

	FILE *fp = fopen(str, "w");

	fputs("Title = \"Instantaneous profile\"\n", fp);
	fputs("variables = \"y\", \"U\", \"V\", \"W\", \"P\", \"R11\", \"R22\", \"R33\", \"R12\", \"R23\", \"R13\", \"Rpu\", \"Rpv\", \"Rpw\", \"Rpp\", \"NU\"\n", fp);
	fprintf(fp, "zone t = \"%i\", i = %i\n", tstep, ms.Ny+1);

	for (int j=0; j<=ms.Ny; j++) {
		fprintf(fp, "%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",
			ms.yc(j),
			um_[j], vm_[j], wm_[j], pm_[j],
			r11_[j], r22_[j], r33_[j],
			r12_[j], r23_[j], r13_[j],
			rpu_[j], rpv_[j], rpw_[j],
			rpp_[j],
			num_[j]	);
	}

	fclose(fp);
}


void Statis::WriteLogfile(const char *path, int tstep, double time, const double mpg[3]) const
{
	FILE *fp;
	char str[1024];
	long int pos = GetLogpos(path, tstep);

	sprintf(str, "%sLOG.dat", path?path:"");

	if (pos > 0) {
		fp = fopen(str, "r+");
		fseek(fp, 0, SEEK_END);
		if (pos < ftell(fp)) {
			char *buf = new char [pos / sizeof(char)];
			fseek(fp, 0, SEEK_SET);
			fread(buf, pos, 1, fp);
			fclose(fp);
			fp = fopen(str, "w");
			fwrite(buf, pos, 1, fp);
			delete[] buf;
		}
	}
	else {
		fp = fopen(str, "w");
		fputs("Title = \"Running log\"\n", fp);
		fputs("variables = \"n\", \"t\", \"ener\", \"taub21\", \"taub22\", \"taub23\", \"mpg1\", \"mpg2\", \"mpg3\", \"Um\", \"Vm\", \"Wm\", \"div\", \"divposx\", \"divposy\", \"divposz\", \"cfl\", \"cflposx\", \"cflposy\", \"cflposz\"\n", fp);
		fputs("zone t = \"statis\"\n", fp);
	}

	fprintf(fp, "%i\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%i\t%i\t%i\t%.18e\t%i\t%i\t%i\n",
		tstep, time, ener_,
		taub_[0], taub_[1], taub_[2],
		mpg  [0], mpg  [1], mpg  [2],
		velm_[0], velm_[1], velm_[2],
		div_, divpos_[0], divpos_[1], divpos_[2],
		cfl_, cflpos_[0], cflpos_[1], cflpos_[2]	);

	fclose(fp);
}


long int Statis::GetLogpos(const char *path, int tstep)
/* find the first line whose time step >= tstep, return the beginning position and record the time of this line */
{
	FILE *fp;
	char str[1024];
	long int pos = 0;
	double n = 0;
	double t = 0.0;

	sprintf(str, "%sLOG.dat", path?path:"");

	if ( ( fp = fopen(str, "r") ) ) {
		fgets(str, 1024, fp);	// skip header
		fgets(str, 1024, fp);	// skip header
		fgets(str, 1024, fp);	// skip header

		while ( n < tstep ) {
			pos = ftell(fp);
			if ( ! fgets(str, 1024, fp) ) break; // reach the end of file
			sscanf(str, "%lf\t%lf", &n, &t);	// if str is empty line or spaces, n and t will not be assigned
		}
		fclose(fp);
	}
	
	return pos;
}

double Statis::GetLogTime(const char *path, int tstep, double mpg[3])
{
	long int pos = GetLogpos(path, tstep);

	char str[1024];
	sprintf(str, "%sLOG.dat", path?path:"");

	double n;
	double t = 0;
	double ener, tau21, tau22, tau23;

	if (pos > 0) {
		FILE *fp = fopen(str, "r");
		fseek(fp, pos, SEEK_SET);
		if(fgets(str, 1024, fp))
			sscanf(str, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				&n, &t, &ener, &tau21, &tau22, &tau23, mpg, mpg+1, mpg+2);
		fclose(fp);
	}

	return t;
}









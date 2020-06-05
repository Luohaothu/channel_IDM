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

	const Scla &u = vel[1], &nuz = vis.SeeVec(3);
	const Scla &v = vel[2], &nuc = vis.SeeScl();
	const Scla &w = vel[3], &nux = vis.SeeVec(1);

	double tau21 = 0;
	double tau22 = 0;
	double tau23 = 0;

	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

		// mean velocity on real boundary
		double um = .5/ms.hy(j) * (u.meanxz(j) * ms.dy(j-1) + u.meanxz(ms.jma(j)) * ms.dy(j));
		double wm = .5/ms.hy(j) * (w.meanxz(j) * ms.dy(j-1) + w.meanxz(ms.jma(j)) * ms.dy(j));
		double vm = v.meanxz(j);

		for (int k=1; k<ms.Nz; k++) {
		for (int i=1; i<ms.Nx; i++) {

			int id =              ms.idx(i,j,k);
			int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
			int im, jm, km;       ms.imx(i,j,k,im,jm,km);
			double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
			double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
			double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);

			double s21 = .5 * ((u[id]-u[jm]) / hyc + (v[id]-v[im]) / hxc);
			double s23 = .5 * ((w[id]-w[jm]) / hyc + (v[id]-v[km]) / hzc);
			double s22 = j==1? (v[jp]-v[id]) / dyc : (v[id]-v[jm]) / dym;

			double r21 = (.5/hyc * (u[id]*dym + u[jm]*dyc) - um)
			           * (.5/hxc * (v[id]*dxm + v[im]*dxc) - vm);
			double r23 = (.5/hyc * (w[id]*dym + w[jm]*dyc) - wm)
			           * (.5/hzc * (v[id]*dzm + v[km]*dzc) - vm);
			double r22 = pow(v[id] - vm, 2.);

			tau21 += (j==1 ? 1 : -1) * (2.*nuz[id] * s21 - r21);
			tau23 += (j==1 ? 1 : -1) * (2.*nux[id] * s23 - r23);
			tau22 += 1./hyc * (nuc[id]*dym + nuc[jm]*dyc) * s22 - r22;
		}}
	}

	taub[0] = .5/((ms.Nx-1)*(ms.Nz-1)) * tau21;
	taub[1] = .5/((ms.Nx-1)*(ms.Nz-1)) * tau22;
	taub[2] = .5/((ms.Nx-1)*(ms.Nz-1)) * tau23;
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
	
	p = fld.SeeScl();
	p.SetLyr(p.SeeLyr(1), 0);
	p.SetLyr(p.SeeLyr(ms.Ny-1), ms.Ny);

	nuc = vis.SeeScl();

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
	int n = 0;
	double t = 0.0;

	sprintf(str, "%sLOG.dat", path?path:"");

	if ( ( fp = fopen(str, "r") ) ) {
		fgets(str, 1024, fp);	// skip header
		fgets(str, 1024, fp);	// skip header
		fgets(str, 1024, fp);	// skip header

		while ( n < tstep ) {
			pos = ftell(fp);
			if ( ! fgets(str, 1024, fp) ) break; // reach the end of file
			sscanf(str, "%i\t%lf", &n, &t);	// if str is empty line or spaces, n and t will not be assigned
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

	int n;
	double t = 0;
	double ener, tau21, tau22, tau23;

	if (pos > 0) {
		FILE *fp = fopen(str, "r");
		fseek(fp, pos, SEEK_SET);
		if(fgets(str, 1024, fp))
			sscanf(str, "%i\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				&n, &t, &ener, &tau21, &tau22, &tau23, mpg, mpg+1, mpg+2);
		fclose(fp);
	}

	return t;
}









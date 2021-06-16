#include <fstream>
#include "Solver.h"
#include "Interp.h"
#include "Filter.h"
#include "IDM.h"
#include "SGS.h"
#include "WM.h"
#include "OFW.h"
#include "Statis.h"

using namespace std;


Solver::Solver(const char *path, bool ifinit):
para(path),
geom(para.Nx, para.Ny, para.Nz, para.Lx, para.Ly, para.Lz),
ms  (geom),
bc  (ms),
sbc (ms),
fld (ms),
fldh(ms),
vis (ms),
fb  (ms),
mpg (3,0.),
time(0.),
step(0)
{
	// initiate mesh grid coordinates
	if (para.prd == 010) geom = Geometry_prdxz(geom).Init(para.dy_min);
	if (para.prd == 110) geom = Geometry_prdz (geom).Init(para.dy_min);

	if (! ifinit) return;

	// initiate the mean pressure gradient with default values, can be adjusted in Evolve
	if (para.prd == 010) { mpg[0] = -1.; mpg[1] = mpg[2] = 0; }
	if (para.prd == 110) { mpg[0] = mpg[1] = mpg[2] = 0; }

	// initiate the flow field for computation
	if (para.inread) {
		ContinueCase();
	} else {
		if (para.dy_min >= 0) InitFieldChan(bc, sbc, para.inener);
		if (para.dy_min <  0) InitFieldBlyr(bc, sbc, para.inener);
	}

	ShowInfo();

	// write initial field if not in-place continuing
	if (step == 0) Output();
}


// ***** initiation ***** //

void Solver::ContinueCase()
{
	if (! para.isInplaceContinue()) {
		// interpolate from another case
		Solver slv_(para.inpath, false);
		slv_.para.inpath[0] = 0;
		slv_.para.inread = para.inread;
		slv_.ContinueCase();
		InitFieldFrom(slv_.fld);
	} else {
		// in-place continue from the current case
		fld.ReadField(para.fieldpath, para.inread, "");
		step = para.inread;
		time = Statis::GetLog(para.statpath, step, mpg);
	}
}

void Solver::InitFieldChan(Boundaries &bc, Boundaries &sbc, double energy)
// initiate flow field (U, V, W, P, including all boundaries) from laminar with random fluctions
{
	Scla &u = fld.GetVec(1);
	Scla &v = fld.GetVec(2);
	Scla &w = fld.GetVec(3);

	fld.InitRand(energy);

	// impose parabolic profile from laminar flow
	for (int j=0; j<=ms.Ny; j++)
		u.AddLyr(1.5 * ms.yc(j) * (2-ms.yc(j)), j);

	// modify flow rate
	u *= 1./u.MeanU(); // bulk mean U = 1.0 due to non-dimensionalization
	w +=  - w.MeanW(); // note: boundaries are modified through using bulk functions

	Bcond::ChannelNoSlip(bc, sbc, ms);
	SetBoundaries(fld.GetVec(), bc, sbc);
}

void Solver::InitFieldBlyr(Boundaries &bc, Boundaries &sbc, double energy)
// initiate flow field (U, V, W, P, including all boundaries) from laminar with random fluctions
{
	Scla &u = fld.GetVec(1);
	Scla &v = fld.GetVec(2);
	Scla &w = fld.GetVec(3);

	fld.InitRand(energy);

	double Ufree = 1.;

	// impose initial profile
	for (int j=0; j<=ms.Ny; j++) {

		double um = Ufree / (PI/2) * atan(100*ms.yc(j)/ms.Ly);

		// damp fluctuations in free stream
		u.MltLyr(1-um/Ufree, j); if (j>0)
		v.MltLyr(1-um/Ufree, j);
		w.MltLyr(1-um/Ufree, j);

		u.AddLyr(um, j);
	}

	Bcond::TblDevelop(bc, sbc, fld.SeeVec(), fld.SeeVec(), Ufree, 1);	
	SetBoundaries(fld.GetVec(), bc, sbc);
}

void Solver::InitFieldFrom(const Flow &fld0)
// initiate flow field (U, V, W, P, including all boundaries) from interpolation of existing fields
{
	Scla &u = fld.GetVec(1);
	Scla &v = fld.GetVec(2);
	Scla &w = fld.GetVec(3);
	Scla &p = fld.GetScl();

	Interp::InterpBulkU(u, fld0.SeeVec(1));
	Interp::InterpBulkV(v, fld0.SeeVec(2));
	Interp::InterpBulkW(w, fld0.SeeVec(3));
	Interp::InterpBulkA(p, fld0.SeeScl());

	for (int j=0; j<=ms.Ny; j++) {
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {
		u(0,j,k) = p(0,j,k) = 0;
		v(i,0,k) = p(i,0,k) = 0;
		w(i,j,0) = p(i,j,0) = 0;
	}}}
}

// ***** output ***** //

void Solver::ShowInfo() const
{
	para.showPara();

	if (para.inread) {
		if (! para.isInplaceContinue()) {
			cout << endl << "Flow fields initiated from existing fields with parameters:" << endl;
			Para(para.inpath).showPara();
		} else
			cout << endl << "Continue from step " << step << ", time " << time << endl;
	} else
		cout << endl << "Flow fields initiated from laminar." << endl;
}

void Solver::Output()
{	
	if (step == 0) {
		geom.WriteMesh (para.statpath);
		geom.WriteMeshY(para.statpath);
	}
	if (step % para.nwrite == 0) {
		fld .WriteField(para.fieldpath, step, "");
		(fldh *= 1./para.dt).WriteField(para.fieldpath, step, "T");
		fld.WriteTecplot(para.fieldpath, 0, time);
		cout << "Files successfully written for step " << step << endl;
	}
	if (step % para.nprint == 0) {

		Statis stas(ms);

		stas.Check(fld, vis, para.Re, para.dt);
		stas.WriteProfile(para.statpath);
		stas.WriteLogfile(para.statpath, step, time, mpg);
	}
}

// ***** computation ***** //

void Solver::CalcVis(Flow &vis, const Vctr &vel, double Re, int sgstyp)
{
	Scla &nuc = vis.GetScl();

	if (sgstyp <= 1) {
		nuc.Set(1./Re);
	} else {
		static SGS sgs(nuc.ms);

		switch (sgstyp) {
		case 2: sgs.Smargorinsky (nuc, vel, Re, .18); break;
		case 3: sgs.DynamicSmarg (nuc, vel         ); break;
		case 4: sgs.DynamicVreman(nuc, vel, Re     ); break;
		}

		nuc += 1./Re;
	}
	#pragma omp parallel
	vis.CellCenter2Edges();
}


void Solver::CalcFb(Vctr &fb, const vector<double> &mpg)
{
	fb[1].Set(- mpg[0]);
	fb[2].Set(- mpg[1]);
	fb[3].Set(- mpg[2]);
}
void Solver::CalcFb(Vctr &fb, const vector<double> &mpg, const Vctr &f)
{
	CalcFb(fb, mpg);
	fb[1] += f[1];
	fb[2] += f[2];
	fb[3] += f[3];
}
void Solver::CalcFb(Vctr &fb, const vector<double> &mpg, const char *filename)
{
	const Mesh &ms = fb.ms;

	fb[1].Set(- mpg[0]);
	fb[2].Set(- mpg[1]);
	fb[3].Set(- mpg[2]);

	// read distributed force
	vector<double> y, f;

	string line;
	ifstream fin; fin.open(filename);

	while (getline(fin, line)) {

		double y_, f_;
		sscanf(line.c_str(), "%le%le",  &y_, &f_);

		y.push_back(y_);
		f.push_back(f_);
	}

	fin.close();

	// Add force
	#pragma omp parallel for
	for (int j=0; j<=ms.Ny; j++) {

		double yc = ms.yc(j);

		int j0 = Interp::BiSearch(yc, y.data(), 0, y.size()-1);
		int j1 = j0 + 1;

		double a = (y[j1] - yc) / (y[j1] - y[j0]);
		double b = (yc - y[j0]) / (y[j1] - y[j0]);

		fb[1].MltLyr(a * f[j0] + b * f[j1], j); // transformed from inner-scale to outer-scale by -mpg
	}
}


void Solver::CalcMpg(vector<double> &mpg, Vctr &vel, double dt)
{
	// double du = vel[1].meanxz(Interp::BiSearch(1., vel.ms.yc(), 0, vel.ms.Ny+1)) - .99;
	double du = vel[1].MeanU() - 1.;
	double dw = vel[3].MeanW();

	vector<double> mpgref(3,0);
	// solve mpg increment with streamwise flowrate 2.0 and spanwise 0
	mpgref[0] = mpg[0] + du / dt;
	mpgref[2] = mpg[2] + dw / dt;
	CalcMpg(mpg, vel, dt, mpgref);
}
void Solver::CalcMpg(vector<double> &mpg, Vctr &vel, double dt, const vector<double> &mpgref)
// solve the increment of mean pressure gradient at n+1/2 step
// given reference MFR at n+1 step or MPG at n+1/2 step
{
	const Mesh &ms = vel.ms;
	Scla &u = vel[1];
	Scla &w = vel[3];

	double dmpg1 = mpgref[0] - mpg[0];
	double dmpg3 = mpgref[2] - mpg[2];

	double limiter = 100;
	dmpg1 = limiter / (PI/2) * atan(dmpg1 / limiter * (PI/2));
	dmpg3 = limiter / (PI/2) * atan(dmpg3 / limiter * (PI/2));

	// update mpg
	mpg[0] += dmpg1;
	mpg[2] += dmpg3;

	// complement mpg increment that was not included in the velocity update step
	#pragma omp parallel for collapse(3)
	for (int j=1; j<ms.Ny; j++) {
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		if (i>0) u(i,j,k) -= dt * dmpg1;
		if (k>0) w(i,j,k) -= dt * dmpg3;
	}}}

	// note: since BC has been applied before, this function should maintain boundary relations
}

void Solver::CalcMpg(vector<double> &mpg, Vctr &vel, const Vctr &fb, double dt)
{
	const Mesh &ms = vel.ms;
	Scla &u = vel[1];
	Scla &w = vel[3];

	int j = Interp::BiSearch(1., ms.yc(), 0, ms.Ny+1);

	double du = u.meanz(1,j) - 1.; // u.meanxz(j) - 1.;
	double dw = 0;                 // w.MeanW();

	double dmpg1 = du/dt / fb[1](1,j,1) * (mpg[0] ? -mpg[0] : 1.);
	double dmpg3 = dw/dt;

	double limiter = 100;
	dmpg1 = limiter / (PI/2) * atan(dmpg1 / limiter * (PI/2));
	dmpg3 = limiter / (PI/2) * atan(dmpg3 / limiter * (PI/2));

	FILE* fp = fopen("dmpg1.dat", "a");
	fprintf(fp, "%.6e\n", dmpg1);
	fclose(fp);

	// complement mpg increment that was not included in the velocity update step
	#pragma omp parallel for collapse(3)
	for (int j=1; j<ms.Ny; j++) {
	for (int k=0; k<=ms.Nz; k++) {
	for (int i=0; i<=ms.Nx; i++) {

		if (i>0) u(i,j,k) -= dt * dmpg1 * fb[1](1,j,1) / (mpg[0] ? -mpg[0] : 1.); // add a fraction of Fx profile
		if (k>0) w(i,j,k) -= dt * dmpg3;
	}}}

	// update mpg
	mpg[0] += dmpg1;
	mpg[2] += dmpg3;
}


void Solver::SetBoundaries(Vctr &vel, const Boundaries &bc, const Boundaries &sbc)
{
	const Mesh &ms = vel.ms;
	// set non-periodic boundaries first
	if (ms.x(0)) Bcond::SetBoundaryX(vel, bc, sbc);
	// if (ms.z(0)) Bcond::SetBoundaryZ(vel, bc, sbc); // not defined yet
	if (ms.y(0)) Bcond::SetBoundaryY(vel, bc, sbc);
	// set periodic boundaries second
	if (!ms.x(0)) Bcond::SetBoundaryX(vel);
	if (!ms.z(0)) Bcond::SetBoundaryZ(vel);
	// if (!ms.y(0)) Bcond::SetBoundaryY(vel);
}


double Solver::InnerScale() const
{
	double tauw = fabs(mpg[0]) ? fabs(mpg[0]) : Statis::WallStress(fld.SeeVec(), vis);
	double utau = sqrt(tauw);
	double Ret  = utau * para.Re;
	return Ret;
}


void Solver::RollBack(Flow &fld, const Flow &fldh, double dt)
{
	Scla &u = fld.GetVec(1);
	Scla &v = fld.GetVec(2);
	Scla &w = fld.GetVec(3);
	Scla &p = fld.GetScl();

	const Scla &uh = fldh.SeeVec(1);
	const Scla &vh = fldh.SeeVec(2);
	const Scla &wh = fldh.SeeVec(3);
	const Scla &dp = fldh.SeeScl();

	( (u *= 1./dt) -= uh ) *= dt;
	( (v *= 1./dt) -= vh ) *= dt;
	( (w *= 1./dt) -= wh ) *= dt;
	( (p *= 1./dt) -= dp ) *= dt;
}


void Solver::RemoveSpanMean(Vctr &vel)
{
	const Mesh &ms = vel.ms;

	Scla &u = vel[1];
	Scla &v = vel[2];

#pragma omp parallel
{
	double um, *usm = new double[ms.Nx+1];
	double vm, *vsm = new double[ms.Nx+1];

	#pragma omp for
	for (int j=1; j<ms.Ny; j++) {
		um = 0;
		for (int i=1; i<ms.Nx; i++)
			um += (usm[i]=u.MeanAz(i,j)) / (ms.Nx-1);

		usm[ms.Nx] = u.MeanAz(ms.Nz,j);

		for (int k=0; k<=ms.Nz; k++)
		for (int i=1; i<=ms.Nx; i++)
			u(i,j,k) -= usm[i] - um;
	}

	#pragma omp for
	for (int j=2; j<ms.Ny; j++) {
		vm = 0;
		for (int i=1; i<ms.Nx; i++)
			vm += (vsm[i]=v.MeanAz(i,j)) / (ms.Nx-1);

		vsm[0]     = v.MeanAz(0,j);
		vsm[ms.Nx] = v.MeanAz(ms.Nx,j);

		for (int k=0; k<=ms.Nz; k++)
		for (int i=0; i<=ms.Nx; i++)
			v(i,j,k) -= vsm[i] - vm;
	}

	delete[] usm;
	delete[] vsm;
}
}


void Solver::SwapBulk(Vctr &vel)
{
	const Mesh &ms = vel.ms;

	Scla &u = vel[1];
	Scla &v = vel[2];
	Scla &w = vel[3];

	int j1 = 8, j2 = 16;
	int k1 = 8, k2 = 16;
	int i1 = 8, i2 = 16;

	char swap = 'x';

	for (int j=0; j<=ms.Ny; j++) {
	for (int i=0; i<=ms.Nx; i++) {
	for (int k=0; k<=ms.Nz; k++) {

		switch (swap) {
		case 'x':
			if (0 < i && i < i1) {
				u(ms.Nx-i2+i,j,k) = u(i,j,k);
				v(ms.Nx-i2+i,j,k) = v(i,j,k);
				w(ms.Nx-i2+i,j,k) = w(i,j,k);
			} else if (i1 <= i && i < i2) {
				u(i,j,k) = u(ms.Nx-i2+i,j,k);
				v(i,j,k) = v(ms.Nx-i2+i,j,k);
				w(i,j,k) = w(ms.Nx-i2+i,j,k);
			}
		break;
		case 'z':
			if (0 < k && k < k1) {
				u(i,j,ms.Nz-k2+k) = u(i,j,k);
				v(i,j,ms.Nz-k2+k) = v(i,j,k);
				w(i,j,ms.Nz-k2+k) = w(i,j,k);
			} else if (k1 <= k && k < k2) {
				u(i,j,k) = u(i,j,ms.Nz-k2+k);
				v(i,j,k) = v(i,j,ms.Nz-k2+k);
				w(i,j,k) = w(i,j,ms.Nz-k2+k);
			}
		break;
		case 'y':
			if (0 < j && j < j1) {
				u(i,j,k) =  u(i,ms.Ny-j,  k);
				w(i,j,k) =  w(i,ms.Ny-j,  k);
				v(i,j,k) = -v(i,ms.Ny-j+1,k);
			} else if (j1 <= j && j < j2) {
				u(i,ms.Ny-j,  k) =  u(i,j,k);
				w(i,ms.Ny-j,  k) =  w(i,j,k);
				v(i,ms.Ny-j+1,k) = -v(i,j,k);
			}
		break;
		}
	}}}


	// // eliminate divergence by solving spanwise Poisson Eq.
	// ////// TODO: not correct
	// Scla p(ms); p.Set(0.);

	// for (int j=1; j<j1; j++) {

	// 	for (int k=1; k<ms.Nz; k++) {
	// 	for (int i=1; i<ms.Nx; i++) {
	// 		p(i,j,k) = vel.Divergence(i,j,k);
	// 	}}

	// 	p.debug_AsciiOutput("", "p1", j, j+1);

	// 	p.fftz(j,j);

	// 	for (int k=0; k<ms.Nzc; k++)
	// 		cout << ms.kz2(k) << endl;

	// 	for (int k=0; k<ms.Nzc; k++) {
	// 	for (int i=1; i<ms.Nx;  i++) {
	// 		p(i,j,2*k)   *= -1. / ms.kz2(k);
	// 		p(i,j,2*k+1) *= -1. / ms.kz2(k);
	// 	}}

	// 	for (int i=1; i<ms.Nx; i++) {
	// 		p(i,j,0) = 0;
	// 		p(i,j,1) = 0;
	// 	}

	// 	p.ifftz(j,j);

	// 	p.debug_AsciiOutput("", "p2", j, j+1);
	// 	exit(0);

	// 	for (int k=1; k<ms.Nz; k++) {
	// 	for (int i=1; i<ms.Nx; i++) {
	// 		w(i,j,k) -= 1./ms.hz(k) * (p(i,j,k) - p(i,j,ms.kma(k)));
	// 	}}

	// 	// if (j==7) cout << vel.Divergence(5,j,5) << endl;
	// }

	// Bcond::SetBoundaryX(vel);
	// Bcond::SetBoundaryZ(vel);



	// // eliminate divergence by solving planar Poisson Eq.
	// Scla p(ms); p.Set(0.);

	// for (int j=1; j<j1; j++) {

	// 	for (int k=1; k<ms.Nz; k++) {
	// 	for (int i=1; i<ms.Nx; i++) {
	// 		p(i,j,k) = vel.Divergence(i,j,k);
	// 	}}

	// 	p.fftxz(j,j);

	// 	for (int k=0; k<ms.Nz-1; k++) {
	// 	for (int i=0; i<ms.Nxc;  i++) {
	// 		p[ms.idfr(i,j,k)] *= -1. / (ms.kx2(i) + ms.kz2(k));
	// 		p[ms.idfi(i,j,k)] *= -1. / (ms.kx2(i) + ms.kz2(k));
	// 	}}

	// 	p[ms.idfr(0,j,0)] = 0;
	// 	p[ms.idfi(0,j,0)] = 0;

	// 	p.ifftxz(j,j);

	// 	for (int k=1; k<ms.Nz; k++) {
	// 	for (int i=1; i<ms.Nx; i++) {
	// 		u(i,j,k) -= 1./ms.hx(i) * (p(i,j,k) - p(ms.ima(i),j,k));
	// 		w(i,j,k) -= 1./ms.hz(k) * (p(i,j,k) - p(i,j,ms.kma(k)));
	// 	}}
	// }

	// Bcond::SetBoundaryX(vel);
	// Bcond::SetBoundaryZ(vel);


	// // correct Reynolds stress
	// for (int j=1; j<=j1; j++) {

	// 	double rs0= Statis::ReynoldsShearStresses(ms.Ny-j+1, vel)[0];
	// 	double ns0= Statis::MeanVisShearStresses (ms.Ny-j+1, vel, vis)[0];
	// 	double rs = Statis::ReynoldsShearStresses(j, vel)[0];
	// 	double ns = Statis::MeanVisShearStresses (j, vel, vis)[0];
		
	// 	if (j == j1)
	// 		cout << ns << ' ' << ns0 << ' ' << rs << ' ' << rs0 << ' ' << vis.GetVec(3)(5,j,5) << endl;

	// 	// for (int k=0; k<=ms.Nz; k++) {
	// 	// for (int i=1; i<=ms.Nx; i++) {
	// 	// 	vis.GetVec(3)(i,j,k) = fmin(fmax(.5 * (rs+rs0) / vel.ShearStrain(i,j,k)[0], -.1), .1);
	// 	// }}

	// 	double rsclvis = fmin(fmax((rs+rs0)/ns, -50.), 50.);

	// 	vis.GetVec(3).MltLyr(rsclvis + 1, j);

	// 	rs = Statis::ReynoldsShearStresses(j, vel)[0];
	// 	ns = Statis::MeanVisShearStresses (j, vel, vis)[0];

	// 	if (j == j1)
	// 		cout << ns << ' ' << ns0 << ' ' << rs << ' ' << rs0 << ' ' << vis.GetVec(3)(5,j,5) << endl << endl;
	// }
	

	// // sponge
	// u.MltLyr((double)(j-j1)/(j2-j), j).AddLyr(u.SeeLyr(ms.Ny-j), j).MltLyr((double)(j2-j)/(j2-j1), j);
	// w.MltLyr((double)(j-j1)/(j2-j), j).AddLyr(w.SeeLyr(ms.Ny-j), j).MltLyr((double)(j2-j)/(j2-j1), j);
	// v.MltLyr((double)(j-j1)/(j2-j), j).MnsLyr(v.SeeLyr(ms.Ny-j+1), j).MltLyr((double)(j2-j)/(j2-j1), j);
}


void Solver::debug_Output(const char path[]) const
{
	int Ny = ms.Ny;

	// char path[1024] = "debug/";
	{
		char names[4][32] = {"U", "V", "W", "P"};
		fld.SeeVec(1).debug_AsciiOutput(path, names[0], 0, Ny+1);
		fld.SeeVec(2).debug_AsciiOutput(path, names[1], 1, Ny+1);
		fld.SeeVec(3).debug_AsciiOutput(path, names[2], 0, Ny+1);
		fld.SeeScl( ).debug_AsciiOutput(path, names[3], 0, Ny+1);
	}
	{
		char names[4][32] = {"UH", "VH", "WH", "DP"};
		fldh.SeeVec(1).debug_AsciiOutput(path, names[0], 0, Ny+1);
		fldh.SeeVec(2).debug_AsciiOutput(path, names[1], 1, Ny+1);
		fldh.SeeVec(3).debug_AsciiOutput(path, names[2], 0, Ny+1);
		fldh.SeeScl( ).debug_AsciiOutput(path, names[3], 0, Ny+1);
	}
	{
		char names[4][32] = {"NUX", "NUY", "NUZ", "NUC"};
		vis.SeeVec(1).debug_AsciiOutput(path, names[0], 1, Ny+1);
		vis.SeeVec(2).debug_AsciiOutput(path, names[1], 0, Ny+1);
		vis.SeeVec(3).debug_AsciiOutput(path, names[2], 1, Ny+1);
		vis.SeeScl( ).debug_AsciiOutput(path, names[3], 0, Ny+1);
	}
	{
		char names[3][32] = {"FBX", "FBY", "FBZ"};
		fb[1].debug_AsciiOutput(path, names[0], 0, Ny+1);
		fb[2].debug_AsciiOutput(path, names[1], 1, Ny+1);
		fb[3].debug_AsciiOutput(path, names[2], 0, Ny+1);
	}
}




#ifdef TIME_TEST_SLV_

#include <sys/time.h>
int cnt = 0;
struct timeval time0, time1;
double t_vis=0, t_bc=0, t_fb=0, t_idm=0, t_mpg=0, t_setbc=0;

void Solver::evolve(double Re, double dt, int bftype)
{
	set_time(time_+dt);

	cout << ++ cnt << endl;

	gettimeofday(&time0, NULL);	CalcVis(vis_, get_vel(), Re, bftype);		gettimeofday(&time1, NULL);	t_vis += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	gettimeofday(&time0, NULL);	Bcond::ChannelNoSlip(bc_, sbc_, ms);		gettimeofday(&time1, NULL);	t_bc  += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	gettimeofday(&time0, NULL);	CalcFb(fb_, mpg_);							gettimeofday(&time1, NULL);	t_fb  += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	gettimeofday(&time0, NULL);	IDM::calc(fld_,fldh_,vis_,fb_,bc_,sbc_,dt);	gettimeofday(&time1, NULL);	t_idm += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	gettimeofday(&time0, NULL);	CalcMpg(mpg_, get_vel(), get_velh(), dt);	gettimeofday(&time1, NULL);	t_mpg += (time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	gettimeofday(&time0, NULL);	Bcond::SetBoundaryY(get_vel(), bc_, sbc_);
								Bcond::SetBoundaryX(get_vel());
								Bcond::SetBoundaryZ(get_vel());				gettimeofday(&time1, NULL);	t_setbc+=(time1.tv_sec - time0.tv_sec) + 1e-6 * (time1.tv_usec - time0.tv_usec);
	if (cnt == 10) {
		cout << endl << endl
			 << "total:\t"  << 1./cnt * (t_vis+t_bc+t_fb+t_idm+t_mpg+t_setbc) << endl
			 << "t_vis:\t"  << t_vis  / t_idm << endl
			 << "t_bc:\t"   << t_bc   / t_idm << endl
			 << "t_fb:\t"   << t_fb   / t_idm << endl
			 << "t_mpg:\t"  << t_mpg  / t_idm << endl
			 << "t_setbc:\t"<< t_setbc/ t_idm << endl << endl;
		exit(0);
	}
}

#endif





// void Solver::ModifyBoundaryVis(Flow &vis, const Vctr &vel, double tau12)
// {
// 	const Mesh &ms = vis.ms;

// 	Scla &nuz = vis.GetVec(3);

// 	const Scla &u = vel[1];
// 	const Scla &v = vel[2];

// 	for (int j=1; j<=ms.Ny; j+= ms.Ny-1) {
// 	for (int k=1; k<ms.Nz; k++) {
// 	for (int i=1; i<ms.Nx; i++) {

// 		double s12 = .5 * (
// 			1./ms.hx(i) * (v(i,j,k) - v(ms.ima(i),j,k)) +
// 			1./ms.hy(j) * (u(i,j,k) - u(i,ms.jma(j),k)) );

// 		nuz(i,j,k) = fmax(.5*tau12/s12, 0);
// 	}}}
// }

// void Solver::ModifyBoundaryVis(Flow &vis, const Vctr &vel, const Vctr &veldns, double Re)
// {
// 	const Mesh &ms = vis.ms, &ms0 = veldns.ms;

// 	const Scla &u = vel[1], &u0 = veldns[1];
// 	const Scla &v = vel[2], &v0 = veldns[2];
// 	const Scla &w = vel[3], &w0 = veldns[3];

// 	// viscosity & sgs-stress on edges aimed for
// 	Vctr shearsgs(ms);
// 	Scla &nux = vis.GetVec(1);
// 	Scla &nuz = vis.GetVec(3);
// 	const Scla &tau23sgs = shearsgs[2];
// 	const Scla &tau12sgs = shearsgs[1];

// 	// sgs shear stress filtered from resolved velocity field
// 	SGS::SubGridShearStress(shearsgs, veldns);

// 	// for top & bottom real boundaries
// 	for (int j=1; j<=ms.Ny; j+=ms.Ny-1) {

// 		double r12dns = 0; // DNS Reynolds stress
// 		double r12les = 0; // LES resolved Reynolds
// 		double v12sgs = 0; // mean sgs-stress filtered from DNS

// 		// the wall normal position on which stresses act
// 		const double y = ms.y(j);

// 		// ***** rescale the DNS-filtered sgs stress by Reynolds stress defect ***** //

// 		// calculate reference (DNS) Reynolds shear stress
// 		int j3u = Interp::BiSearch(y, ms0.yc(), 0, ms0.Ny), j4u = j3u+1;
// 		int j3v = Interp::BiSearch(y, ms0.y(),  1, ms0.Ny), j4v = j3v+1;

// 		double um3 = u0.meanxz(j3u), y3u = ms0.yc(j3u);
// 		double um4 = u0.meanxz(j4u), y4u = ms0.yc(j4u);
// 		double vm3 = v0.meanxz(j3v), y3v = ms0.y(j3v);
// 		double vm4 = v0.meanxz(j4v), y4v = ms0.y(j4v);

// 		double um = (um3 * (y4u-y) + um4 * (y-y3u)) / (y4u-y3u);
// 		double vm = (vm3 * (y4v-y) + vm4 * (y-y3v)) / (y4v-y3v);

// 		#pragma omp parallel for reduction(+: r12dns)
// 		for (int k=1; k<ms0.Nz; k++) { double z = ms0.zc(k);
// 		for (int i=1; i<ms0.Nx; i++) { double x = ms0.x(i);
// 			// calculate uv on z-edges and interpolate to the desired y position
// 			r12dns += (Interp::InterpNodeU(x,y,z,u0) - um)
// 			        * (Interp::InterpNodeV(x,y,z,v0) - vm)
// 			        * ((ms.Nx-1) * (ms.Nz-1)) / ((ms0.Nx-1) * (ms0.Nz-1)); // rescale to match the other two
// 		}}

// 		// calculate resolved Reynolds stress and mean sgs-stress
// 		um3 = u.meanxz(j3u = ms.jma(j)); y3u = ms.yc(j3u);
// 		um4 = u.meanxz(j4u = j);         y4u = ms.yc(j4u);
// 		vm  = v.meanxz(j);
// 		um  = (um3 * (y4u-y) + um4 * (y-y3u)) / (y4u-y3u);

// 		for (int k=1; k<ms.Nz; k++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hyc=ms.hy(j);
// 		for (int i=1; i<ms.Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxc=ms.hx(i);

// 			int id =        ms.idx(i,j,k);
// 			int im, jm, km; ms.imx(i,j,k,im,jm,km);

// 			v12sgs += tau12sgs[id];
// 			r12les += (.5/hyc * (u[id]*dym + u[jm]*dyc) - um)
// 			        * (.5/hxc * (v[id]*dxm + v[im]*dxc) - vm);
// 		}}

// 		const double rescale1 = fabs((r12dns - r12les) / v12sgs);
// 		const double rescale2 = fabs((y4u-y3u) / (1-fabs(y-1)) / log((1-fabs(y4u-1))/(1-fabs(y3u-1))));

// 		// // rescale physical viscosity by replacing low-order differencing with log law
// 		// double dyU2 = (um4-um3) / log((1-fabs(y4u-1))/(1-fabs(y3u-1))) / (1-fabs(y-1));
// 		// double dyU1 = (um4-um3) / (y4u-y3u);
// 		// const double rescale2 = fabs(dyU2 / dyU1);

// 		// ***** rescale viscosity to account for low-order differencing error ***** //

// 		// construct low- & high-order differenced U gradient
// 		// int j5u = j==1 ? j4u+1 : j3u-1;
// 		// double um5 = u.meanxz(j5u);
// 		// double y5u = ms.yc(j5u);

// 		// double dyU2 = (y4u+y5u-2*y) / (y4u-y3u) / (y3u-y5u) * um3
// 		//             + (y3u+y5u-2*y) / (y3u-y4u) / (y4u-y5u) * um4
// 		//             + (y3u+y4u-2*y) / (y3u-y5u) / (y5u-y4u) * um5;

// 		cout << rescale1 << " " << rescale2 << endl;

// 		// ***** modify physical & eddy viscosity accordingly ***** //

// 		for (int k=1; k<=ms.Nz; k++) {
// 		for (int i=1; i<=ms.Nx; i++) {
// 			// boundary strain rate of the coarse velocity field
// 			const double *sr = vel.ShearStrain(i,j,k);

// 			double s12 = sr[0], tau12 = tau12sgs(i,j,k) * rescale1;
// 			double s23 = sr[1], tau23 = tau23sgs(i,j,k) * rescale1;

// 			if (k<ms.Nz) nuz(i,j,k) = rescale2/Re + fmax(.5 * tau12 / s12, 0);
// 			if (i<ms.Nx) nux(i,j,k) = rescale2/Re + fmax(.5 * tau23 / s23, 0);
// 		}}
// 	}
// }






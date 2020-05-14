#pragma once

#include "Basic.h"


class Boundary
{
public:
	Boundary(const Mesh &ms):
	ms(ms),
	Nx(ms.Nx),
	Ny(ms.Ny),
	Nz(ms.Nz)
	{
		b1_ = new double[(Ny+1) * (Nz+1)]; // back
		b2_ = new double[(Ny+1) * (Nz+1)]; // front
		b3_ = new double[(Nx+1) * (Nz+1)]; // bottom
		b4_ = new double[(Nx+1) * (Nz+1)]; // top
		b5_ = new double[(Nx+1) * (Ny+1)]; // left
		b6_ = new double[(Nx+1) * (Ny+1)]; // right
	};

	~Boundary()
	{
		delete[] b1_; delete[] b2_;
		delete[] b3_; delete[] b4_;
		delete[] b5_; delete[] b6_;
	};

	double& b1(int j, int k) { return b1_[j*(Nz+1) + k]; };
	double& b2(int j, int k) { return b2_[j*(Nz+1) + k]; };
	double& b3(int i, int k) { return b3_[k*(Nx+1) + i]; };
	double& b4(int i, int k) { return b4_[k*(Nx+1) + i]; };
	double& b5(int i, int j) { return b5_[j*(Nx+1) + i]; };
	double& b6(int i, int j) { return b6_[j*(Nx+1) + i]; };

	double b1(int j, int k) const { return b1_[j*(Nz+1) + k]; };
	double b2(int j, int k) const { return b2_[j*(Nz+1) + k]; };
	double b3(int i, int k) const { return b3_[k*(Nx+1) + i]; };
	double b4(int i, int k) const { return b4_[k*(Nx+1) + i]; };
	double b5(int i, int j) const { return b5_[j*(Nx+1) + i]; };
	double b6(int i, int j) const { return b6_[j*(Nx+1) + i]; };

private:
	const Mesh &ms;
	const int Nx, Ny, Nz;
	double *b1_, *b2_;
	double *b3_, *b4_;
	double *b5_, *b6_;

};


class Boundaries
{
public:
	Boundaries(const Mesh &ms):
	ub_(ms),
	vb_(ms),
	wb_(ms)
	{};

	double& ub1(int j, int k) { return ub_.b1(j,k); };
	double& ub2(int j, int k) { return ub_.b2(j,k); };
	double& ub3(int i, int k) { return ub_.b3(i,k); };
	double& ub4(int i, int k) { return ub_.b4(i,k); };
	double& ub5(int i, int j) { return ub_.b5(i,j); };
	double& ub6(int i, int j) { return ub_.b6(i,j); };

	double& vb1(int j, int k) { return vb_.b1(j,k); };
	double& vb2(int j, int k) { return vb_.b2(j,k); };
	double& vb3(int i, int k) { return vb_.b3(i,k); };
	double& vb4(int i, int k) { return vb_.b4(i,k); };
	double& vb5(int i, int j) { return vb_.b5(i,j); };
	double& vb6(int i, int j) { return vb_.b6(i,j); };

	double& wb1(int j, int k) { return wb_.b1(j,k); };
	double& wb2(int j, int k) { return wb_.b2(j,k); };
	double& wb3(int i, int k) { return wb_.b3(i,k); };
	double& wb4(int i, int k) { return wb_.b4(i,k); };
	double& wb5(int i, int j) { return wb_.b5(i,j); };
	double& wb6(int i, int j) { return wb_.b6(i,j); };

	double ub1(int j, int k) const { return ub_.b1(j,k); };
	double ub2(int j, int k) const { return ub_.b2(j,k); };
	double ub3(int i, int k) const { return ub_.b3(i,k); };
	double ub4(int i, int k) const { return ub_.b4(i,k); };
	double ub5(int i, int j) const { return ub_.b5(i,j); };
	double ub6(int i, int j) const { return ub_.b6(i,j); };

	double vb1(int j, int k) const { return vb_.b1(j,k); };
	double vb2(int j, int k) const { return vb_.b2(j,k); };
	double vb3(int i, int k) const { return vb_.b3(i,k); };
	double vb4(int i, int k) const { return vb_.b4(i,k); };
	double vb5(int i, int j) const { return vb_.b5(i,j); };
	double vb6(int i, int j) const { return vb_.b6(i,j); };

	double wb1(int j, int k) const { return wb_.b1(j,k); };
	double wb2(int j, int k) const { return wb_.b2(j,k); };
	double wb3(int i, int k) const { return wb_.b3(i,k); };
	double wb4(int i, int k) const { return wb_.b4(i,k); };
	double wb5(int i, int j) const { return wb_.b5(i,j); };
	double wb6(int i, int j) const { return wb_.b6(i,j); };

private:
	Boundary ub_;
	Boundary vb_;
	Boundary wb_;
};


namespace Bcond
{

//apply boundary conditions for velocity fields
void SetBoundaryX(Vctr &vel, const Boundaries &bc, const Boundaries &sbc);
void SetBoundaryY(Vctr &vel, const Boundaries &bc, const Boundaries &sbc);
void SetBoundaryX(Vctr &vel);
void SetBoundaryZ(Vctr &vel);

// set homogeneous boundary for cell-centered scalars
void SetBoundaryX(Scla &q, int ord);
void SetBoundaryY(Scla &q, int ord);
void SetBoundaryZ(Scla &q, int ord);

// Channel: no-slip on both (real) walls
void ChannelNoSlip(Boundaries &bc, Boundaries &sbc, const Mesh &ms);

// Channel off-wall: Dirichlet on both (virtual) boundaries
void ChannelDirichlet(Boundaries &bc, Boundaries &sbc, const Mesh &ms, const Vctr &vel);

// Channel: homogeneous Robin BC ( u - l du/dy = 0 ) on both walls
void ChannelRobin(Boundaries &bc, Boundaries &sbc, const Mesh &ms, const Vctr &vel);

// TBL: no-slip on bottom wall; u=U, dv/dy=dw/dy=0 on top boundary
void TblCycling(Boundaries &bc, Boundaries &sbc, const Mesh &ms, double ufree);

} // namespace Bcond












// void SetBoundaryY(Vctr &vel, Vctr &velh,
// 	const Boundaries &bc, const Boundaries &sbc, const Mesh &ms)
// {
// 	int Nx = ms.Nx;
// 	int Ny = ms.Ny;
// 	int Nz = ms.Nz;

// 	Scla &u = vel(1), &uh = velh(1);
// 	Scla &v = vel(2), &vh = velh(2);
// 	Scla &w = vel(3), &wh = velh(3);

// 	for (int k=0; k<=Nz; k++) {
// 	for (int i=0; i<=Nx; i++) {

// 		double ubc3 = bc.ub3(i,k), sbu3 = sbc.ub3(i,k);
// 		double vbc3 = bc.vb3(i,k), sbv3 = sbc.vb3(i,k);
// 		double wbc3 = bc.wb3(i,k), sbw3 = sbc.wb3(i,k);

// 		double ubc4 = bc.ub4(i,k), sbu4 = sbc.ub4(i,k);
// 		double vbc4 = bc.vb4(i,k), sbv4 = sbc.vb4(i,k);
// 		double wbc4 = bc.wb4(i,k), sbw4 = sbc.wb4(i,k);

// 		uh(i,0,k) = 1./dt * (ubc3 - sbu3 * u(i,1,k) - u(i,0,k));
// 		vh(i,1,k) = 1./dt * (vbc3 - sbv3 * v(i,2,k) - v(i,1,k));
// 		wh(i,0,k) = 1./dt * (wbc3 - sbw3 * w(i,1,k) - w(i,0,k));

// 		uh(i,Ny,k) = 1./dt * (ubc4 - sbu4 * u(i,Ny-1,k) - u(i,Ny,k));
// 		vh(i,Ny,k) = 1./dt * (vbc4 - sbv4 * v(i,Ny-1,k) - v(i,Ny,k));
// 		wh(i,Ny,k) = 1./dt * (wbc4 - sbw4 * w(i,Ny-1,k) - w(i,Ny,k));

// 		u(i,0,k) += dt * uh(i,0,k);
// 		v(i,1,k) += dt * vh(i,1,k);
// 		w(i,0,k) += dt * wh(i,0,k);

// 		u(i,Ny,k) += dt * uh(i,Ny,k);
// 		v(i,Ny,k) += dt * vh(i,Ny,k);
// 		w(i,Ny,k) += dt * wh(i,Ny,k);
// 	}}
// };








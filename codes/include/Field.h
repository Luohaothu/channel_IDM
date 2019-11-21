# pragma once


# include "Basic.h"
# include "IDM.h"
# include "SGS.h"


# include <iostream>
using namespace std;


class Field: private Mesh
{
	public:
		Field(const Mesh &mesh):
			Mesh(mesh), idm(mesh), sgs(mesh),
			U(mesh), UH(mesh), UBC(Mesh(Nx,2-1,Nz,Lx,Ly,Lz)),
			P(mesh), DP(mesh, true), NU(mesh),
			mpg{0,0,0}
		{
			// UBC.meshGet().y[0] = UBC.meshGet().yc[0] = y[1];
			// UBC.meshGet().y[1] = UBC.meshGet().yc[1] = y[Ny];
			UBC.meshGet().y[0] = 0;
			UBC.meshGet().y[1] = 0;
			UBC.meshGet().yc[0] = y[1];
			UBC.meshGet().yc[1] = y[Ny];
		};

		Vctr U, UH, UBC;// UH also serve to store the RHS of delta u^* (increment of intermediate velocities), and delta u^* itself
		Scla P, DP, NU;	// dp also serve to store the RHS of the pressure Poisson equation
		double mpg[3];

		// initiation
		void initField(double energy);	// initiate flow fields from laminar with random fluctuations
		void initField(Field &field);	// initiate flow fields from existing fields

		// boundary (UBC) process
		void bcond(int tstep);
		void bcond(Vctr &U0);

		// computation related
		void getnu(double Re, int bftype);
		void getup(double dt, int nthrds);
		void applyBC();
		void applyBC(double dt);
		void removeSpanMean();

		// IO functions
		void writeField(char *path, int tstep);
		void readField (char *path, int tstep);
		void writeTecplot(char *path, int tstep, double time);
		void debug_Output(int tstep);

	private:
		IDM idm;
		SGS sgs;
};





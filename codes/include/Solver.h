# pragma once

# include "Basic.h"
# include "IDM.h"
# include "SGS.h"
# include "DA.h"


class Solver
{
	private:
		const Mesh &mesh;
		IDM idm;
		SGS sgs;
		int nthrds;

	public:
		Feld FLD, FLDH, VIS, BC;
		Vctr FB; double mpg[3];


	public:
		Solver(const Mesh &mesh, const Mesh &bmesh):
			mesh(mesh), idm(mesh), sgs(mesh),
			FLD(mesh), FLDH(mesh), VIS(mesh), BC(bmesh),
			FB(mesh), mpg{0,0,0} {};

		void config(int n);
		void debug_Output(int tstep);

		/***** default cumputation *****/
		void evolve(double Re, double dt, int bftype)
		{
			getbc();
			getnu(Re, bftype);
			getfb();
			getup(dt);
			fixfr(dt);
			if (bftype == 1) removeSpanMean(); // for MFU
		};
		/***** off-wall BC computation *****/
		void evolve_ofw(double Re, double dt, int bftype, const Vctr &U, const double MPG[3])
		{
			getbc(U);
			getnu(Re, bftype);
			getfb();
			getup(dt);
			fixfr(dt, MPG);
		};
		/***** data assimilation *****/
		void assimilate(DA &da, double dt)
		{
			int n = 10; double e = 1e-6, a = .01; // DA parameters
			// int n = 10000; double e = 1e-4, a = .1; // DA parameters

			// // store the unassimilated flow field for recovery later
			// Feld FLD_temp(mesh), FLDH_temp(mesh);
			// FLD_temp.V[1] = FLD.V[1]; FLDH_temp.V[1] = FLDH.V[1];
			// FLD_temp.V[2] = FLD.V[2]; FLDH_temp.V[2] = FLDH.V[2];
			// FLD_temp.V[3] = FLD.V[3]; FLDH_temp.V[3] = FLDH.V[3];
			// FLD_temp.S    = FLD.S;    FLDH_temp.S    = FLDH.S;


			while (da.ifIter(FLD.V, e, n)) {     // converged or reached max iterations yet
				rollback(dt);              // roll back the flow fields to the old time step
				da.getAdj(FLD.V, VIS, dt); // compute adjoint state using the new time step fields
				getfb(da.getForce(a));     // apply the assimilating force
				getup(dt);                 // solve the time step again under the new force
			}

			// // recover the flow field to the unassimilated state
			// FLD.V[1] = FLD_temp.V[1]; FLDH.V[1] = FLDH_temp.V[1];
			// FLD.V[2] = FLD_temp.V[2]; FLDH.V[2] = FLDH_temp.V[2];
			// FLD.V[3] = FLD_temp.V[3]; FLDH.V[3] = FLDH_temp.V[3];
			// FLD.S    = FLD_temp.S;    FLDH.S    = FLDH_temp.S;
		};

	private:
		// construct boundary conditions
		void getbc();
		void getbc(const Vctr &U);
		// construct viscosity
		void getnu(double Re, int bftype);
		// construct body forces
		void getfb();
		void getfb(const Vctr &F);
		// evolve velocity and pressure fields by 1 time step
		void getup(double dt);
		// modify flow rates by adjusting mean pressure gradients
		void fixfr(double dt);
		void fixfr(double dt, const double MPG[3]);
		// data manipulation
		void rollback(double dt);
		void removeSpanMean();

};




		




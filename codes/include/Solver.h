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

		Solver(const Mesh &mesh, const Mesh &bmesh):
			mesh(mesh), idm(mesh), sgs(mesh),
			FLD(mesh), FLDH(mesh), VIS(mesh), BC(bmesh),
			FB(mesh), mpg{0,0,0} {};

		void config(int n);
		void debug_Output(int tstep);


	/***** default cumputation *****/
	public:
		void evolve(double Re, double dt, int bftype)
		{
			getbc();
			getnu(Re, bftype);
			getfb();
			getup(dt);
			if (bftype == 1) removeSpanMean(); // for MFU
		};
	private:
		void getbc() { BC.reset(); };		// construct boundary conditions
		void getnu(double Re, int bftype);	// construct viscosity
		void getfb();						// construct body forces
		void getup(double dt)				// evolve velocity and pressure fields by 1 time step
			{ idm.calc(FLD, FLDH, mpg, VIS, FB, BC, dt); };
		void removeSpanMean();				// data manipulation


	/***** off-wall BC computation *****/
	public:
		void evolve_ofw(double Re, double dt, int bftype, const Vctr &U, const double MPG[3])
		{
			getbc(U);
			getnu(Re, bftype);
			getfb();
			getup(dt);

			// keep the OFW mean pressure gradients equal to FC
			double dmpg1 = MPG[0] - mpg[0];
			double dmpg3 = MPG[2] - mpg[2];
			mpg[0] += dmpg1;
			mpg[2] += dmpg3;
			for (int j=1; j<mesh.Ny; j++) {
				FLD.V[1].lyrAdd(- dt * dmpg1, j);
				FLD.V[3].lyrAdd(- dt * dmpg3, j);
			}
		};
	private:
		void getbc(const Vctr &U);


	/***** data assimilation computation *****/
	public:
		void assimilate(double dt, DA &da)
		{
			int n = 10; double e = 1e-4, a = .1; // DA parameters

			while (da.ifIter(FLD.V, e, n)) {     // converged or reached max iterations yet
				da.getAdj(FLD.V, VIS, dt); // compute adjoint state using the new time step fields
				getfb(da.getForce(a));     // apply the assimilating force
				rollback(dt);              // roll back the flow fields to the old time step
				getup(dt);                 // solve the time step again under the new force
			}
		};
	private:
		void getfb(const Vctr &F);
		void rollback(double dt);

};




		




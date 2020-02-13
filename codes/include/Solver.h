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
		void evolve_ofw(double Re, double dt, int bftype, const Vctr &U)
		{
			getbc(U);
			getnu(Re, bftype);
			getfb();
			getup(dt);
		};
	private:
		void getbc(const Vctr &U);


	/***** data assimilation computation *****/
	public:
		void evolve_da(double Re, double dt, int bftype, DA &da)
		{
			getbc();
			getnu(Re, bftype);
			getfb();
			getup(dt);

			while (da.ifIter(FLD.V, 0.01, 10)) {
				// compute adjoint state using the new time step fields
				da.calcAdj(FLD.V, VIS, dt);
				// apply the assimilating force
				getfb(da.getForce(0.1));
				// roll back the flow fields to the old time step
				FLD.V[1] -= (FLDH.V[1] *= dt);
				FLD.V[2] -= (FLDH.V[2] *= dt);
				FLD.V[3] -= (FLDH.V[3] *= dt);
				FLD.S    -= (FLDH.S    *= dt);
				// solve the time step again under the new force
				getup(dt);
			}
		};
	private:
		void getfb(const Vctr &F) {
			(FB[1] = F[1]) += -mpg[0];
			(FB[2] = F[2]) += -mpg[1];
			(FB[3] = F[3]) += -mpg[2];
		};

};




		




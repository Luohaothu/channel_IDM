# pragma once

# include "Basic.h"
# include "IDM.h"
# include "SGS.h"
# include "Interp.h"


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

		// void initiate(double energy) { FLD.initrand(energy); };
		// void initiate(const Feld &f) { FLD.initfrom(f); };

		void config(int n)
		{
			int tempn = nthrds;
			if (tempn != (nthrds = idm.ompset(n)))
				std::cout << "\nnumber of OpenMP threads: " << nthrds << '\n' << std::endl;
		};

		void evolve(double Re, double dt, int bftype)
		{
			getbc();
			getnu(Re, bftype);
			getfb();
			getup(dt);
			if (bftype == 1) removeSpanMean(); // for MFU
		};

		// construct boundary conditions
		void getbc() { BC.reset(); };
		void getbc(const Vctr &U)
		{
			const Mesh &ms0 = U.meshGet();
			int jb = (ms0.Ny - mesh.Ny) / 2;

			Interp(U[1], BC.V[1]).layerPrdLin(jb,        0);
			Interp(U[1], BC.V[1]).layerPrdLin(ms0.Ny-jb, 1);

			Interp(U[2], BC.V[2]).layerPrdLin(jb+1,      0);
			Interp(U[2], BC.V[2]).layerPrdLin(ms0.Ny-jb, 1);

			Interp(U[3], BC.V[3]).layerPrdLin(jb,        0);
			Interp(U[3], BC.V[3]).layerPrdLin(ms0.Ny-jb, 1);

			// Interp(U[1], BC.V[1]).bulkFilter('U');
			// Interp(U[2], BC.V[2]).bulkFilter('X');
			// Interp(U[3], BC.V[3]).bulkFilter('U');
		};

		// construct viscosty
		void getnu(double Re, int bftype)
		{
			Scla &NU = VIS.S; const Vctr &U = FLD.V;

			switch (bftype) {
			case 2: sgs.smargorinsky (NU, U, Re, .18); NU += 1./Re; break;
			case 3: sgs.dynamicsmarg (NU, U);          NU += 1./Re; break;
			case 4: sgs.dynamicvreman(NU, U, Re);      NU += 1./Re; break;
			default:                                   NU  = 1./Re;
			}

			VIS.CC2EG();
		};

		// construct body forces
		void getfb()
		{
			FB[1] = - mpg[0];
			FB[2] = - mpg[1];
			FB[3] = - mpg[2];
		};

		// evolve velocity and pressure fields by 1 time step
		void getup(double dt)
		{
			idm.calc(FLD, FLDH, mpg, VIS, FB, BC, dt);
		};

		// data manipulation
		void removeSpanMean();
		

		void debug_Output(int tstep);

};


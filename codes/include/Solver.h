# pragma once

# include "Field.h"
# include "IDM.h"
# include "SGS.h"


class Solver
{
	public:
		Solver(Mesh &mesh, Field &field):
		mesh(mesh),
		field(field),
		idm(mesh, field),
		sgs(mesh, field.NU, field.U)
		{};

		void config(int n)
		{
			int tempn = nthrds;
			if (tempn != (nthrds = idm.ompset(n)))
				std::cout << "\nnumber of OpenMP threads: " << nthrds << std::endl << std::endl;
		};

		void getup(double Re, double dt, int bftype)
		/* evolve velocity and pressure fields by 1 time step */
		{
			this->getbc();
			this->getnu(Re, bftype);
			this->getfb();
			this->idm.upcalc(dt);
			if (bftype == 1) this->removeSpanMean(); // for MFU
		};

		void getup(double Re, double dt, int bftype, Vctr &U0)
		/* evolve velocity and pressure fields by 1 time step */
		{
			this->getbc(U0);
			this->getnu(Re, bftype);
			this->getfb();
			this->idm.upcalc(dt);
			if (bftype == 1) this->removeSpanMean(); // for MFU
		};

	private:
		Mesh &mesh;
		Field &field;

		IDM idm;
		SGS sgs;

		int nthrds;

		void getbc();          // construct boundary conditions
		void getbc(Vctr &U0);
		void getnu(double Re, int bftype); // construct viscosty
		void getfb();          // construct body forces
		void removeSpanMean(); // data manipulation
};
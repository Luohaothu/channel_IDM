# pragma once

# include "Basic.h"

class Press
{
	public:

		void U2P(Scla &P, Scla &NU, Vctr &U, Vctr &UT, double dt)
		{
			Vctr V(U.meshGet()); // here V is the old time step and U is the new time step, different to that in rhs()
			V.com1.bulkCpy(UT.com1).bulkMlt(-dt);
			V.com2.bulkCpy(UT.com2).bulkMlt(-dt);
			V.com3.bulkCpy(UT.com3).bulkMlt(-dt);

			for (int j=0; j<=U.meshGet().Ny; j++) {
				V.com1.layersAdd(U.com1, j, j);
				V.com2.layersAdd(U.com2, j, j);
				V.com3.layersAdd(U.com3, j, j);
			}

			this->rhs(P, NU, V, U, dt);
			this->poisson(P);
		};
		void rhs(Scla &R, Scla &NU, Vctr &U, Vctr &V, double dt);
		void poisson(Scla &P);
};
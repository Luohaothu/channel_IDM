# pragma once

# include "Basic.h"

class Press
{
	public:
		static void u2p(Feld FLD, const Vctr &UT, const Feld &VIS, double dt)
		{
			Feld RHS(FLD.meshGet());
			rhs(RHS, FLD.V, UT, VIS, dt);
			poisson(FLD.S = RHS.S);
		};

		static void rhs(Feld &R, const Vctr &U, const Vctr &UT, const Feld &VIS, double dt);
		static void poisson(Scla &P);
};












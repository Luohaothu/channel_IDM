# pragma once

# include "Basic.h"

class SGS: private Mesh
{
	public:
		SGS(const Mesh &mesh): Mesh(mesh) {};
		void smargorinsky(Scla &EV, Vctr &U, double Cs, double Re);
		void dynamicsmarg(Scla &EV, Vctr &U);
		void dynamicvreman(Scla &EV, Vctr &U, double Re);
};
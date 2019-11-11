# pragma once

# include "Basic.h"

class SGS: private Mesh
{
	public:
		SGS(const Mesh &mesh): Mesh(mesh) {};

	// 	void evcalc(Scla &EV, Vctr &U, double Re=0, double Cs=0);
		
	// private:
	// 	double Cs; // smargorinsky constant
		void smargorinsky(Scla &EV, Vctr &U, double Cs, double Re);
		void dynamicsmarg(Scla &EV, Vctr &U);
};
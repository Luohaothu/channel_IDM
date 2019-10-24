# pragma once

# include "Basic.h"

class SGS: private Mesh
{
	public:
		SGS(const Mesh &mesh): Mesh(mesh) {};

		void evcalc(Scla &EV, Vctr &U, double Re, double Cs);
		
	private:
		double Cs; // smargorinsky constant
		void smargorinsky(Scla &NU, Vctr &U, double Re);
};
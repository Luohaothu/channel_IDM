# pragma once

# include "Basic.h"

class SGS: private Mesh
{
	public:
		SGS(const Mesh &mesh, Scla &EV, Vctr &U):
		Mesh(mesh), EV(EV), U(U) {};
		
		void smargorinsky(double Cs, double Re);
		void dynamicsmarg();
		void dynamicvreman(double Re);

	private:
		Scla &EV;
		Vctr &U;
};
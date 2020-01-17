# pragma once

# include "Basic.h"

class SGS: private Mesh
{
	public:
		SGS(const Mesh &mesh): Mesh(mesh) {};
		
		void smargorinsky (Scla &EV, const Vctr &U, double Re, double Cs);
		void dynamicsmarg (Scla &EV, const Vctr &U);
		void dynamicvreman(Scla &EV, const Vctr &U, double Re);
};
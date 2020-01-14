# pragma once

# include "Basic.h"

class Field: private Mesh
{
	public:
		Field(const Mesh &mesh);
		~Field() {};

		Vctr U, UH, FB, UBC; // UH also serve to store the RHS of delta u^* (increment of intermediate velocities), and delta u^* itself
		Scla P, DP, NU, PBC; // dp also serve to store the RHS of the pressure Poisson equation
		double mpg[3];

		// initiation
		void initField(double energy); // initiate flow fields from laminar with random fluctuations
		void initField(Field &field);  // initiate flow fields from existing fields

		// IO functions
		void writeField(char *path, int tstep);
		Field& readField(char *path, int tstep);
		void writeTecplot(char *path, int tstep, double time);
		void debug_Output(int tstep);

	private:
		void reset();   // set irrelavent variables (all but U,V[1:],W,P) to 0
};





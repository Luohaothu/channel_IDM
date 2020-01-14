# pragma once

# include "Basic.h"
# include "Field.h"

class IDM: private Mesh
{
	public:
		IDM(const Mesh &mesh, Field &field);
		// IDM(const Mesh &mesh,
		// 	Vctr &U, Vctr &UH, Vctr &FB, Vctr &UBC,
		// 	Scla &P, Scla &DP, Scla &NU, Scla &PBC, double (&mpg)[3]);

		// configuration
		int ompset(int n);

		// computation interface
		void upcalc(double dt);

	private:
		Vctr &U, &UH, &FB, &UBC;
		Scla &P, &DP, &NU, &PBC;
		double (&mpg)[3];

		Scla NUX, NUY, NUZ;
		double *u, *uh, *ubc, *fbx, *nux,
		       *v, *vh, *vbc, *fby, *nuy,
		       *w, *wh, *wbc, *fbz, *nuz,
		       *p, *dp, *fdp, *nu;

		// computations functions
		       
		// step 1: interpolate viscosity from cell-centers to edges
		void visset();
		// step 2: calculate RHS of momentum equations
		void urhs1();
		void urhs2();
		void urhs3();
		void mbc();
		// step 3: calculate intermedia velocities (solve TDMAs)
		void getuh1(double dt);
		void getuh2(double dt);
		void getuh3(double dt);
		// step 4: calculate projector (perform FFTs)
		void rhsdp(double dt);
		void getfdp(double refp);
		// step 5: update velocity & pressure fields
		void update(double dt);
		void meanpg(double dt);
		// step 6: update doundaries
		void applyBC(double dt);
		void pressBD();          // boundary conditions for pressure
		void veldtBD(double dt); // modify boundary of velocity-time-derivative with given BC
		void velocBD();          // boundary conditions for velocities
};






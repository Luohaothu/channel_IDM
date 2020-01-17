# pragma once

# include <omp.h>
# include "Basic.h"


class IDM: private Mesh
{
	public:
		IDM(const Mesh &mesh): Mesh(mesh) {};

		// configuration
		int ompset(int n);

		// computation interface
		void calc(Feld &FLD, Feld &FLDH, double mpg[3],
			const Feld &VIS, const Vctr &FB, const Feld &BC, double dt)
		{
			uhcalc(FLDH.V, FLD, VIS, BC, FB, dt);
			dpcalc(FLDH, FLD.S, BC.V, dt);
			upcalc(FLD, FLDH, mpg, dt);
			applyBC(FLD, FLDH, BC, dt);
		};

	private:
		// computations functions
		void uhcalc(Vctr &UH,
			const Feld &FLD, const Feld &VIS, const Feld &BC, const Vctr &FB, double dt)
		{
			# pragma omp parallel
			{
				urhs1(UH[1].blkGet(), FLD, VIS, FB[1].blkGet());
				urhs2(UH[2].blkGet(), FLD, VIS, FB[2].blkGet());
				urhs3(UH[3].blkGet(), FLD, VIS, FB[3].blkGet());
				mbc(UH, FLD, VIS, BC);
				getuh1(UH, FLD.V, VIS, dt);
				getuh2(UH, FLD.V, VIS, dt);
				getuh3(UH, FLD.V, VIS, dt);
			}
		};
		void dpcalc(Feld &FLDH,
			const Scla &P, const Vctr &UBC, double dt)
		{
			Scla &DP = FLDH.S;
			rhsdp(DP.blkGet(), FLDH.V, UBC, dt); // rdp (which shares memory with dp)
			DP.fft();                            // rdp->frdp
			getfdp(DP.blkGetF(), P.av(1));       // frdp->fdp
			DP.ifft();                           // fdp->dp
		};
		void upcalc(Feld &FLD, Feld &FLDH, double mpg[3], double dt)
		{
			update(FLD, FLDH, dt);
			meanpg(FLD.V, mpg, dt);
		};
		void applyBC(Feld &FLD, Feld &FLDH, const Feld &BC, double dt)
		{
			// // extrapolate UBC from real boundary to virtual boundary at new time step
			// UBC[1].lyrMlt(2.*h[1]/dy[0], 0).lyrMns(U[1][1], 0).lyrMlt(dy[0]/dy[1], 0);
			// UBC[3].lyrMlt(2.*h[1]/dy[0], 0).lyrMns(U[3][1], 0).lyrMlt(dy[0]/dy[1], 0);
			// UBC[1].lyrMlt(2.*h[Ny]/dy[Ny], 1).lyrMns(U[1][Ny-1], 1).lyrMlt(dy[Ny]/dy[Ny-1], 1);
			// UBC[3].lyrMlt(2.*h[Ny]/dy[Ny], 1).lyrMns(U[3][Ny-1], 1).lyrMlt(dy[Ny]/dy[Ny-1], 1);
			
			// when UBC is aligned to the virtual boundary
			pressBD(FLD.S, FLDH.S, BC.S);
			veldtBD(FLDH.V, FLD.V, BC.V, dt);
			velocBD(FLD.V, BC.V);
		};

		// subroutines for computation
		// step 1: calculate RHS of momentum equations
		void urhs1(double *ruh, const Feld &FLD, const Feld &VIS, double *fbx);
		void urhs2(double *rvh, const Feld &FLD, const Feld &VIS, double *fby);
		void urhs3(double *rwh, const Feld &FLD, const Feld &VIS, double *fbz);
		void mbc  (Vctr &UH, const Feld &FLD, const Feld &VIS, const Feld &BC);
		// step 2: calculate intermedia velocities (solve TDMAs)
		void getuh1(Vctr &UH, const Vctr &U, const Feld &VIS, double dt);
		void getuh2(Vctr &UH, const Vctr &U, const Feld &VIS, double dt);
		void getuh3(Vctr &UH, const Vctr &U, const Feld &VIS, double dt);
		// step 3: calculate projector (perform FFTs)
		void rhsdp(double *rdp, const Vctr &UH, const Vctr &UBC, double dt);
		void getfdp(double *fdp, double refp);
		// step 4: update velocity & pressure fields
		void update(Feld &FLD, Feld &FLDH, double dt);
		void meanpg(Vctr &U, double mpg[3], double dt);
		// step 5: update doundaries
		void pressBD(Scla &P, Scla &DP, const Scla &PBC);                  // boundary conditions for pressure
		void veldtBD(Vctr &UH, const Vctr &U, const Vctr &UBC, double dt); // modify boundary of velocity-time-derivative with given BC
		void velocBD(Vctr &U, const Vctr &UBC);                            // boundary conditions for velocities
};






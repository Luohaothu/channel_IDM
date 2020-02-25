# pragma once

# include "Basic.h"


class DA
{
	public:
		DA(const Mesh &mesh): _F(mesh), _FLDH(mesh), _UE(mesh), _MSK(mesh)
		{
			_iter = _erro = 0; _vexp = setMask(mesh);
			_F.reset(); _FLDH.reset(); _UE.reset();
		};

		bool getExp(double time, const Vctr &UE);
		bool ifIter(const Vctr &U, double e, int n);
		void getAdj(const Vctr &U, const Feld &VIS, double dt);
		const Vctr& getForce(double alpha);

	private:
		Vctr _F;        // driving force to assimilate the flow
		Feld _FLDH;     // adjoint velocities & pressure
		Vctr _UE, _MSK; // experiment velocity field & the mask function
		int _iter;      // iteration number
		double _erro;   // residual of the iteration
		double _vexp;   // volume of the space where experiment data exits

		int Nx, Ny, Nz, Nxz;
		double dx, dz, dx2, dz2, *dy, *h;
		int *ipa, *ima, *kpa, *kma;

		const Mesh& setMesh(const Mesh &ms)
		{
			Nx = ms.Nx; Ny = ms.Ny; Nz = ms.Nz; Nxz = ms.Nxz;
			dx = ms.dx; dz = ms.dz; dx2 = ms.dx2; dz2 = ms.dz2;
			dy = ms.dy; h = ms.h;
			ipa = ms.ipa; ima = ms.ima; kpa = ms.kpa; kma = ms.kma;
			return ms;
		};

		double setMask(const Mesh &ms)
		{
			_MSK.reset(1.);
			_MSK[1].lyrSet(0., 0).lyrSet(0., ms.Ny);
			_MSK[2].lyrSet(0., 1).lyrSet(0., ms.Ny).lyrSet(0., 0);
			_MSK[3].lyrSet(0., 0).lyrSet(0., ms.Ny);

			// _MSK.reset();
			// for (int j=0; j<2; j++) {
			// 	_MSK[1].lyrSet(1., j+1).lyrSet(1., ms.Ny-(j+1));
			// 	_MSK[2].lyrSet(1., j+2).lyrSet(1., ms.Ny-(j+1));
			// 	_MSK[3].lyrSet(1., j+1).lyrSet(1., ms.Ny-(j+1));
			// }

			for (int j=1; j<ms.Ny; j++) {
			for (int k=0; k<ms.Nz; k++) {
			for (int i=0; i<ms.Nx; i++) {
				_FLDH.S.id(i,j,k) = _MSK.module(i,j,k);
			}}}
			return _FLDH.S.bulkMeanU();
		};

		void urhs  (Vctr &UH, const Vctr &U, const Vctr &UE, const Vctr &MSK);
		void getuh1(Vctr &UH, const Vctr &U, const Feld &VIS, double dt);
		void getuh2(Vctr &UH, const Vctr &U, const Feld &VIS, double dt);
		void getuh3(Vctr &UH, const Vctr &U, const Feld &VIS, double dt);
		void rhsdp (Scla &DP, const Vctr &UH, double dt);
		void getfdp(Scla &DP, double refp);
		void update(Vctr &UH, const Scla &DP, double dt);
};








# pragma once

# include "Basic.h"


class DA
{
	public:
		DA(const Mesh &mesh): _F(mesh), _FLDH(mesh), _UE(mesh), _MSK(mesh)
		{
			_MSK.reset(1.);
			_MSK[1].lyrSet(0., 0).lyrSet(0., mesh.Ny);
			_MSK[2].lyrSet(0., 1).lyrSet(0., mesh.Ny).lyrSet(0., 0);
			_MSK[3].lyrSet(0., 0).lyrSet(0., mesh.Ny);
		};

		bool getExp(double time, const Vctr &UE);
		bool ifIter(const Vctr &U, double e, int n);

		void calcAdj(const Vctr &U, const Feld &VIS, double dt)
		{
			urhs(_FLDH.V, U, _UE, _MSK);

			getuh1(_FLDH.V, U, VIS, dt);
			getuh2(_FLDH.V, U, VIS, dt);
			getuh3(_FLDH.V, U, VIS, dt);

			rhsdp(_FLDH.S, _FLDH.V, dt);
			_FLDH.S.fft();
			getfdp(_FLDH.S, 0.);
			_FLDH.S.ifft();

			update(_FLDH.V, _FLDH.S, dt);

			// homogeneous BC applied automatically
		};

		const Vctr& getForce(double lambda)
		{
			_F[1] += (_FLDH.V[1] *= lambda);
			_F[2] += (_FLDH.V[2] *= lambda);
			_F[3] += (_FLDH.V[3] *= lambda);
			_iter ++;
			return _F;
		};

	private:
		Vctr _F;        // driving force to assimilate the flow
		Feld _FLDH;     // adjoint velocities & pressure
		Vctr _UE, _MSK; // experiment velocity field & the mask function
		int _iter;      // iteration number

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

		void urhs  (Vctr &UH, const Vctr &U, const Vctr &UE, const Vctr &MSK);
		void getuh1(Vctr &UH, const Vctr &U, const Feld &VIS, double dt);
		void getuh2(Vctr &UH, const Vctr &U, const Feld &VIS, double dt);
		void getuh3(Vctr &UH, const Vctr &U, const Feld &VIS, double dt);
		void rhsdp (Scla &DP, const Vctr &UH, double dt);
		void getfdp(Scla &DP, double refp);
		void update(Vctr &UH, const Scla &DP, double dt);
};








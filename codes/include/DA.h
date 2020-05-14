#pragma once

#include "Basic.h"


class DA
{
	public:
		DA(const Mesh &ms);

		bool GetExp(double time, const Vctr &vele);
		void Reset();
		bool IfIter(double en, const Vctr &vel);
		const Vctr& GetAsmForce(const Vctr &vel, const Flow &vis, double dt, double alpha);
		
		void WriteLog(double time);
		void WriteForce(const char *path, int tstep) const;
		
	private:
		Vctr fb_;     // driving force to assimilate the flow
		Flow fldh_;   // adjoint velocities & pressure
		Vctr vele_;   // experiment velocity field
		Vctr msk_;    // mask function indicating positions where experiment data are available
		double erro_; // residual of the iteration
		int iter_;    // iteration number

		static void SetMsk(Vctr &msk);

		static void CalcAdjoint(Flow &fldh, const Vctr &vel, const Flow &vis, double dt);
		static void CalcForce(Vctr &fb, const Vctr &velh, double alpha);

		static double urhs(Flow &fldh, const Vctr &vel, const Vctr &vele, const Vctr &msk);
		static void getuh1(Vctr &velh, const Vctr &vel, const Flow &vis, double dt);
		static void getuh2(Vctr &velh, const Vctr &vel, const Flow &vis, double dt);
		static void getuh3(Vctr &velh, const Vctr &vel, const Flow &vis, double dt);
		static void rhsdp (Scla &rdp, const Vctr &velh, double dt);
		static void getfdp(Scla &fdp, double refp);
		static void update(Vctr &velh, const Scla &dp, double dt);
};








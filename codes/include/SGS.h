#pragma once

#include "Basic.h"


namespace SGS
{

void Smargorinsky (Scla &nut, const Vctr &vel, double Re, double Cs);
void DynamicSmarg (Scla &nut, const Vctr &vel);
void DynamicVreman(Scla &nut, const Vctr &vel, double Re);

void SubGridStress(Vctr &shear, Vctr &normal, const Vctr &vel);
void SubGridShearStress(Vctr &shear, const Vctr &vel);

} // namespace SGS

// note: MUST NOT filter and assign values to an array in the same transverse



// class SGS
// {
// public:
// 	SGS(const Mesh &ms);

// 	void Smargorinsky (Scla &nut, const Vctr &vel, double Re, double Cs);
// 	void DynamicSmarg (Scla &nut, const Vctr &vel);
// 	void DynamicVreman(Scla &nut, const Vctr &vel, double Re);

// private:
// 	Scla s11, s22, s33, s12, s23, s13;
// 	Scla m11, m22, m33, m12, m23, m13;
// 	Scla l11, l22, l33, l12, l23, l13;
// };



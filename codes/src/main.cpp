# include <iostream>

# include "Para.h"
# include "Mesh.h"
# include "Field.h"
# include "Statis.h"
# include "IDM.h"

using namespace std;


void output(int tstep, double time, class Para *ppara, class Mesh *pmesh, class Field *pfield, class Statis *pstat);


int main()
{
	class Para	*ppara	= new class Para	("XINDAT");
	class Mesh	*pmesh	= new class Mesh	(ppara->dim());
	class Field	*pfield	= new class Field	(ppara->dim());
	class Statis*pstat	= new class Statis	(ppara->dim());
	class IDM	*pIDM	= new class IDM		();

	int tstep = 0;
	double time = 0.0;
	
	// initialization
	pmesh->initMesh(ppara->len(), ppara->statpath, ppara->dy_min);
	pfield->initField(ppara->init_ener, pmesh);
	pIDM->initIDM(ppara->Re, ppara->dt, pmesh);
	output(tstep, time, ppara, pmesh, pfield, pstat);

	while (tstep++ < ppara->Nt) {
		time += ppara->dt;

		// set boundary conditions
		pfield->bcond(tstep);

		// time evolution
		pIDM->uhcalc(pfield->UH(), pfield->U(), pfield->P(), pfield->UBC());
		pIDM->dpcalc(pfield->DP(), pfield->UH(), pfield->UBC(), pfield);
		pIDM->upcalc(pfield->U(), pfield->P(), pfield->UH(), pfield->DP(), pmesh, pfield);
		pfield->applyBC();

		// output
		output(tstep, time, ppara, pmesh, pfield, pstat);
	}

}


void output(int tstep, double time, class Para *ppara, class Mesh *pmesh, class Field *pfield, class Statis *pstat)
{
	if (tstep % ppara->nwrite == 0)	{
		pfield->writeField(ppara->fieldpath, tstep);
		pfield->writeTecplot(ppara->fieldpath, tstep, time, pmesh);
		// pfield->debug_Output(tstep);
	}
	if (tstep % ppara->nprint == 0) {
		pstat->checkStat(pfield->UP(), pmesh, pfield, ppara->dt);
		pstat->writeProfile(ppara->statpath, -1, pmesh->yc);
		pstat->writeLogfile(ppara->statpath, tstep, time);
	}
}











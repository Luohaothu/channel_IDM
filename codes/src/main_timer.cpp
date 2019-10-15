# include <iostream>
# include <stdlib.h>
# include <string.h>

# include <sys/time.h>

# include "Para.h"
# include "Mesh.h"
# include "Field.h"
# include "Statis.h"
# include "IDM.h"

using namespace std;


void output(int tstep, double time, class Para *ppara, class Mesh *pmesh, class Field *pfield, class Statis *pstat);
void input(int *tstep, double *time, class Field *pfield, class Para *ppara, class Mesh *pmesh);


int main()
{
	class Para	*ppara	= new class Para	("XINDAT");	ppara->showPara();
	class Mesh	*pmesh	= new class Mesh	(ppara->dim());
	class Field	*pfield	= new class Field	(ppara->dim());
	class Statis*pstat	= new class Statis	(ppara->dim());
	class IDM	*pIDM	= new class IDM		();

	int tstep = 0;	double time = 0.0;

	struct timeval *time0 = new struct timeval;
	struct timeval *time1 = new struct timeval;
	struct timeval *time2 = new struct timeval;
	long t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0;
	
	// initialization
	pmesh->initMesh(ppara->len(), ppara->statpath, ppara->dy_min);
	if (! pmesh->checkYmesh(ppara->statpath, pmesh->y)) exit(0);

	pIDM->initIDM(ppara->Re, ppara->dt, pmesh);

	// pfield->bcond(0);
	if (! ppara->nread) pfield->initField(ppara->init_ener, pmesh);
	else input(& tstep, & time, pfield, ppara, pmesh);

	if (tstep == 0)	output(0, 0.0, ppara, pmesh, pfield, pstat);


gettimeofday(time0, NULL);

	// main loop
	while (tstep++ < ppara->Nt) {
		time += ppara->dt;

		// set boundary conditions
		pfield->bcond(tstep);

		// time evolution

gettimeofday(time1, NULL);
		pIDM->ruhcalc(pfield->UH, pfield->U, pfield->P, pfield->UBC);
gettimeofday(time2, NULL);
t1 += 1000000 * (time2->tv_sec - time1->tv_sec) + time2->tv_usec - time1->tv_usec;

gettimeofday(time1, NULL);
		pIDM->uhcalc(pfield->UH, pfield->U);
gettimeofday(time2, NULL);
t2 += 1000000 * (time2->tv_sec - time1->tv_sec) + time2->tv_usec - time1->tv_usec;

gettimeofday(time1, NULL);
		pIDM->dpcalc(pfield->DP, pfield->UH, pfield->UBC, pfield);
gettimeofday(time2, NULL);
t3 += 1000000 * (time2->tv_sec - time1->tv_sec) + time2->tv_usec - time1->tv_usec;

gettimeofday(time1, NULL);
		pIDM->upcalc(pfield->U, pfield->P, pfield->UPH, pmesh);
gettimeofday(time2, NULL);
t4 += 1000000 * (time2->tv_sec - time1->tv_sec) + time2->tv_usec - time1->tv_usec;
		
gettimeofday(time1, NULL);
		// pfield->applyBC();
		pfield->applyBC(ppara->dt);
gettimeofday(time2, NULL);
t5 += 1000000 * (time2->tv_sec - time1->tv_sec) + time2->tv_usec - time1->tv_usec;

gettimeofday(time1, NULL);
		pfield->bodyForce(ppara->bftype);
gettimeofday(time2, NULL);
t6 += 1000000 * (time2->tv_sec - time1->tv_sec) + time2->tv_usec - time1->tv_usec;


		// output
		output(tstep, time, ppara, pmesh, pfield, pstat);
	}


gettimeofday(time1, NULL);
t7 = time1->tv_sec - time0->tv_sec;

	cout<< "\nuhcalc : " << (double)t2 / t1 << endl
		<< "\ndpcalc : " << (double)t3 / t1 << endl
		<< "\nupcalc : " << (double)t4 / t1 << endl
		<< "\nbccalc : " << (double)t5 / t1 << endl
		<< "\nbfcalc : " << (double)t6 / t1 << endl
		<< "\nstep average : " << (double)t7 / ppara->Nt << " s\n" << endl;

}



void output(int tstep, double time, class Para *ppara, class Mesh *pmesh, class Field *pfield, class Statis *pstat)
{
	if (tstep % ppara->nwrite == 0)	{
		pfield->writeField(ppara->fieldpath, tstep);
		pfield->writeFieldDt(ppara->fieldpath, tstep); // incorrect at tstep 0
		// pfield->writeTecplot(ppara->fieldpath, tstep, time, pmesh);
		// pfield->debug_Output(tstep);
		cout << "Files successfully written for step " << tstep << endl;
	}
	if (tstep % ppara->nprint == 0) {
		pstat->checkStat(pfield->UP, pmesh, pfield, ppara->Re, ppara->dt);
		pstat->writeProfile(ppara->statpath, -1, pmesh->yc);
		pstat->writeLogfile(ppara->statpath, tstep, time);
	}
}


void input(int *tstep, double *time, class Field *pfield, class Para *ppara, class Mesh *pmesh)
{
	char fn1[1024], fn2[1024], fn3[1024];
	sprintf(fn1, "%s%s", ppara->inpath, "XINDAT");

	class Para *ppara0 = new class Para (fn1);
	class Mesh *pmesh0 = new class Mesh (ppara0->dim());
	class Field *pfield0 = new class Field (ppara0->dim());
	class Statis *pstat0 = new class Statis (ppara0->dim());
	sprintf(fn2, "%s%s", ppara->inpath, ppara0->statpath);
	sprintf(fn3, "%s%s", ppara->inpath, ppara0->fieldpath);

	cout << "\nParameters for continue computing:" << endl;
	ppara0->showPara();

	if (! strcmp(ppara->inpath, "")) pstat0->readLogfile(fn2, ppara->nread, time);
	*tstep = *time < 1e-10 ? 0 : ppara->nread; // continue last case or start a new case depending on whether a continue time is read

	pmesh0->initMesh(ppara0->len(), fn2, ppara0->dy_min);
	pfield0->readField(fn3, ppara->nread);
	pfield->initField(pfield0, pmesh0, pmesh);

	delete ppara0;
	delete pmesh0;
	delete pfield0;
	delete pstat0;
}










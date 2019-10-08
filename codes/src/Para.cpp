# include <iostream>
# include <stdio.h>
# include <string.h>

# include "Para.h"


using namespace std;



Para::Para(string filename)
{
	/* if no filename is provided, initiate with default parameters */

	// file paths
	sprintf(fieldpath, "fielddata/");
	sprintf(probepath, "probedata/");
	sprintf(statpath, "statdata/");
	sprintf(postpath, "postdata/");

	// physical parameters
	Re = 100.0;
	init_ener = 0.5;
	bftype = 0;

	// grid settings
	Nx = 4;	// peridoic without overlap
	Ny = 5;	// from wall to wall
	Nz = 4;	// peridoic without overlap
	Lx = 6.2832;
	Ly = 2.0;
	Lz = 3.1416;
	dy_min = 0.2;

	// time evolution control
	Nt = 1000000;
	dt = 1e-3;

	// output control
	nwrite = 10000;
	nprint = 1000;
	nprobe = 0;
	jprbs = new int;
	jprbs[0] = 0;

	// input control
	nread = 0;
	sprintf(inpath, "");

	/* if filename provided, read the file */
	if (filename != "")	this->readPara(filename);
}

void Para::readPara(string filename)
{
	char str[1024], *s;
	FILE *fp = fopen(filename.c_str(), "r");

	while ( fgets(str, 1024, fp) ) {
		if ( (s = strstr(str, "//")) )	{ *s = '\0'; }	// strip the comments

		if ( strstr(str, "fieldpath") )	{ s = strchr(str, '='); s = strchr(++s, '\"'); sscanf(++s, "%[^\"]", fieldpath); }
		if ( strstr(str, "probepath") )	{ s = strchr(str, '='); s = strchr(++s, '\"'); sscanf(++s, "%[^\"]", probepath); }
		if ( strstr(str, "statpath") )	{ s = strchr(str, '='); s = strchr(++s, '\"'); sscanf(++s, "%[^\"]", statpath); }
		if ( strstr(str, "postpath") )	{ s = strchr(str, '='); s = strchr(++s, '\"'); sscanf(++s, "%[^\"]", postpath); }

		if ( strstr(str, "Re") )	{ s = strchr(str, '='); sscanf(++s, "%lf", & Re); }
		if ( strstr(str, "init_ener") )	{ s = strchr(str, '='); sscanf(++s, "%lf", & init_ener); }
		if ( strstr(str, "bftype") ){ s = strchr(str, '='); sscanf(++s, "%i", & bftype); }
		if ( strstr(str, "Nx") )	{ s = strchr(str, '='); sscanf(++s, "%i", & Nx); }
		if ( strstr(str, "Ny") )	{ s = strchr(str, '='); sscanf(++s, "%i", & Ny); }
		if ( strstr(str, "Nz") )	{ s = strchr(str, '='); sscanf(++s, "%i", & Nz); }
		if ( strstr(str, "Lx") )	{ s = strchr(str, '='); sscanf(++s, "%lf", & Lx); }
		if ( strstr(str, "Ly") )	{ s = strchr(str, '='); sscanf(++s, "%lf", & Ly); }
		if ( strstr(str, "Lz") )	{ s = strchr(str, '='); sscanf(++s, "%lf", & Lz); }
		if ( strstr(str, "dy_min") ){ s = strchr(str, '='); sscanf(++s, "%lf", & dy_min); }
		if ( strstr(str, "Nt") )	{ s = strchr(str, '='); sscanf(++s, "%i", & Nt); }
		if ( strstr(str, "dt") )	{ s = strchr(str, '='); sscanf(++s, "%lf", & dt); }
		if ( strstr(str, "nwrite") ){ s = strchr(str, '='); sscanf(++s, "%i", & nwrite); }
		if ( strstr(str, "nprint") ){ s = strchr(str, '='); sscanf(++s, "%i", & nprint); }
		if ( strstr(str, "nprobe") ){ s = strchr(str, '='); sscanf(++s, "%i", & nprobe); }
		if ( strstr(str, "jprbs") )	{ jprbs = this->parseJprbs(str); }
		if ( strstr(str, "nread") )	{ s = strchr(str, '='); sscanf(++s, "%i", & nread); }
		if ( strstr(str, "inpath") ){ s = strchr(str, '='); s = strchr(++s, '\"'); sscanf(++s, "%[^\"]", inpath); }
	}
}

int* Para::parseJprbs(char *str)
{
	char *s;
	int *js, jn, j;

	s = strchr(str, '=');
	sscanf(++s, "%i", & jn);

	js = new int [jn + 1];
	js[0] = jn;

	for (j=1; j<=jn; j++) {
		if ( (s = strchr(s, ',')) )	sscanf(++s, "%i", & js[j]);
		else {
			js[0] = j-1;
			break;
		}
	}

	return js;
}

void Para::showPara()
{
	cout << "\n-------------------- PARAMETERS --------------------" << endl;
	cout << "File paths:" << endl;
	cout << "fieldpath = \"" << fieldpath << '\"' << endl;
	cout << "probepath = \"" << probepath << '\"' << endl;
	cout << "statpath = \"" << statpath << '\"' << endl;
	cout << "postpath = \"" << postpath << '\"' << endl;
	cout << "\nPhysical Parameters:" << endl;
	cout << "Re = " << Re << ", init_ener = " << init_ener << endl;
	cout << "body force type = " << (bftype==1 ? "MFU" : bftype==2 ? "LES" : "DNS") << endl;
	cout << "\nGrid Settings:" << endl;
	cout << "Nx = " << Nx << ", Ny = " << Ny << ", Nz = " << Nz << endl;
	cout << "Lx = " << Lx << ", Ly = " << Ly << ", Lz = " << Lz << endl;
	cout << "dy_min = " << dy_min << endl;
	cout << "\nTime Control:" << endl;
	cout << "Nt = " << Nt << ", dt = " << dt << endl;
	cout << "\nOutput control:" << endl;
	cout << "nwrite = " << nwrite << ", nprint = " << nprint << ", nprobe = " << nprobe << endl;
	cout << "jprbs = ";
		for (int j=1; j<=jprbs[0]; j++)	cout << jprbs[j] << ", ";
		cout << "total " << jprbs[0] << " layers to be probed." << endl;
	cout << "\nInput control:" << endl;
	cout << "nread = " << nread << endl;
	cout << "inpath = \"" << inpath << '\"' << endl;
	cout << "----------------------------------------------------" << endl;
}



// # define DEBUG // g++ -I include src/Para.cpp
# ifdef DEBUG
int main()
{
	class Para *ppara0 = new class Para();
	class Para *ppara = new class Para("XINDAT");
}
# endif

















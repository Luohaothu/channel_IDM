# include <iostream>
# include <stdio.h>
# include <string.h>
# include <sys/stat.h>

# include "Para.h"


using namespace std;


Para::Para(const char *path)
{
	XINDAT_modify_time = 0;
	
	/* initiate default parameters */
	{
		strcpy(workpath, "");

		// file paths
		strcpy(fieldpath, "fielddata/");
		strcpy(probepath, "probedata/");
		strcpy(statpath,  "statdata/" );
		strcpy(postpath,  "postdata/" );
		strcpy(inpath,    ""          );

		// computation control
		bftype = 0;
		nthrds = 0;

		// physical parameters
		Re = 100.0;
		inener = 0.5;

		// grid settings
		Nx = 5;	// peridoic without overlap
		Ny = 5;	// from wall to wall
		Nz = 5;	// peridoic without overlap
		Lx = 6.2832;
		Ly = 2.0;
		Lz = 3.1416;
		dy_min = 0.2;
		prd = 010;

		// time evolution control
		Nt = 1000000;
		dt = 1e-3;

		// IO control
		inread = 0;
		nwrite = 10000;
		nprint = 1000;
		nprobe = 0;
		jprbs[0] = 0;
	}

	/* if workpath provided, read the input file */
	if (path != NULL) {
		strcpy(workpath, path);
		checkPara();
	}
}

void Para::checkPara(int tstep)
{
	char str[1024]; strcat(strcpy(str, workpath), "XINDAT");
	struct stat buf; if (stat(str, &buf)) return; // get file state and store in buf, return if file not exist
	time_t tempt = XINDAT_modify_time;

	if (tempt < (XINDAT_modify_time = buf.st_mtime)) {

		readPara();

		if ((int)tempt != 0) {
			cout << endl << str << " updated at step " << tstep << " as:" << endl;
			showPara();
		}
	}
}


void parseJprbs(int *jprbs, char *str)
{
	char *s = strchr(str, '=');
	sscanf(++s, "%i", & jprbs[0]);
	if (jprbs[0] > 1023) jprbs[0] = 1023;
	for (int j=1; j<=jprbs[0]; j++) {
		if ( (s = strchr(s, ',')) )	sscanf(++s, "%i", & jprbs[j]);
		else { jprbs[0] = j-1; break; }
	}
}

void Para::readPara()
{
	char str[1024], *s;
	FILE *fp = fopen(strcat(strcpy(str, workpath), "XINDAT"), "r");

	while ( fgets(str, 1024, fp) ) {
		if ( (s = strstr(str, "//")) )	{ *s = '\0'; }	// strip the comments

		if (strstr(str, "fieldpath")) { s = strchr(str, '='); s = strchr(++s, '\"'); if (! sscanf(++s, "%[^\"]", fieldpath)) fieldpath[0] = '\0'; }
		if (strstr(str, "probepath")) { s = strchr(str, '='); s = strchr(++s, '\"'); if (! sscanf(++s, "%[^\"]", probepath)) probepath[0] = '\0'; }
		if (strstr(str, "statpath") ) { s = strchr(str, '='); s = strchr(++s, '\"'); if (! sscanf(++s, "%[^\"]", statpath))  statpath [0] = '\0'; }
		if (strstr(str, "postpath") ) { s = strchr(str, '='); s = strchr(++s, '\"'); if (! sscanf(++s, "%[^\"]", postpath))  postpath [0] = '\0'; }
		if (strstr(str, "inpath")   ) { s = strchr(str, '='); s = strchr(++s, '\"'); if (! sscanf(++s, "%[^\"]", inpath))    inpath   [0] = '\0'; }
		if (strstr(str, "bftype")   ) { s = strchr(str, '='); sscanf(++s, "%i",  & bftype); }
		if (strstr(str, "nthrds")   ) { s = strchr(str, '='); sscanf(++s, "%i",  & nthrds); }
		if (strstr(str, "Re")       ) { s = strchr(str, '='); sscanf(++s, "%lf", & Re); }
		if (strstr(str, "inener")   ) { s = strchr(str, '='); sscanf(++s, "%lf", & inener); }
		if (strstr(str, "Nx")       ) { s = strchr(str, '='); sscanf(++s, "%i",  & Nx); }
		if (strstr(str, "Ny")       ) { s = strchr(str, '='); sscanf(++s, "%i",  & Ny); }
		if (strstr(str, "Nz")       ) { s = strchr(str, '='); sscanf(++s, "%i",  & Nz); }
		if (strstr(str, "Lx")       ) { s = strchr(str, '='); sscanf(++s, "%lf", & Lx); }
		if (strstr(str, "Ly")       ) { s = strchr(str, '='); sscanf(++s, "%lf", & Ly); }
		if (strstr(str, "Lz")       ) { s = strchr(str, '='); sscanf(++s, "%lf", & Lz); }
		if (strstr(str, "dy_min")   ) { s = strchr(str, '='); sscanf(++s, "%lf", & dy_min); }
		if (strstr(str, "prd")      ) { s = strchr(str, '='); sscanf(++s, "%i",  & prd); }
		if (strstr(str, "Nt")       ) { s = strchr(str, '='); sscanf(++s, "%i",  & Nt); }
		if (strstr(str, "dt")       ) { s = strchr(str, '='); sscanf(++s, "%lf", & dt); }
		if (strstr(str, "inread")   ) { s = strchr(str, '='); sscanf(++s, "%i",  & inread); }
		if (strstr(str, "nwrite")   ) { s = strchr(str, '='); sscanf(++s, "%i",  & nwrite); }
		if (strstr(str, "nprint")   ) { s = strchr(str, '='); sscanf(++s, "%i",  & nprint); }
		if (strstr(str, "nprobe")   ) { s = strchr(str, '='); sscanf(++s, "%i",  & nprobe); }
		if (strstr(str, "jprbs")    ) { parseJprbs(jprbs, str); }
	}

	fclose(fp);

	// change relative-to-XINDAT paths to be absolute or relative to the executable
	strcpy(fieldpath, strcat(strcpy(str, workpath), fieldpath));
	strcpy(probepath, strcat(strcpy(str, workpath), probepath));
	strcpy(statpath,  strcat(strcpy(str, workpath), statpath));
	strcpy(postpath,  strcat(strcpy(str, workpath), postpath));
	if (inpath[0] != '/' && inpath[0] != '\\')
		strcpy(inpath, strcat(strcpy(str, workpath), inpath));
}


void Para::showPara() const
{
	cout << endl;
	cout << "-------------------- PARAMETERS --------------------" << endl;

	cout << "File paths:" << endl;
	cout << "fieldpath = \"" << fieldpath << '\"' << endl;
	cout << "probepath = \"" << probepath << '\"' << endl;
	cout << "statpath = \"" << statpath << '\"' << endl;
	cout << "postpath = \"" << postpath << '\"' << endl;
	cout << "inpath = \"" << inpath << '\"' << endl;
	cout << endl;
	cout << "Computation control:" << endl;
	cout << "body force type = " << (bftype==0 ? "DNS" : bftype==1 ? "MFU" : bftype==2 ? "LES_SM" : bftype==3 ? "LES_DSM" : "LES_DVR") << endl;
	cout << "number of threads = " << nthrds << endl;
	cout << endl;
	cout << "Physical Parameters:" << endl;
	cout << "Re = " << Re << ", inener = " << inener << endl;
	cout << endl;
	cout << "Grid Settings:" << endl;
	cout << "Nx = " << Nx << ", Ny = " << Ny << ", Nz = " << Nz << endl;
	cout << "Lx = " << Lx << ", Ly = " << Ly << ", Lz = " << Lz << endl;
	cout << "dy_min = " << dy_min << endl;
	cout << "prd = " << prd << endl;
	cout << endl;
	cout << "Time Control:" << endl;
	cout << "Nt = " << Nt << ", dt = " << dt << endl;
	cout << endl;
	cout << "IO control:" << endl;
	cout << "inread = " << inread << endl;
	cout << "nwrite = " << nwrite << endl;
	cout << "nprint = " << nprint << endl;
	cout << "nprobe = " << nprobe << endl;
	cout << "jprbs = ";
		for (int j=1; j<=jprbs[0]; j++)	cout << jprbs[j] << ", ";
		cout << "total " << jprbs[0] << " layers to be probed." << endl;

	cout << "----------------------------------------------------" << endl;
}



// # define DEBUG // g++ -I include src/Para.cpp
#ifdef DEBUG
int main()
{
	class Para *ppara0 = new class Para();
	class Para *ppara = new class Para("XINDAT");
}
#endif

















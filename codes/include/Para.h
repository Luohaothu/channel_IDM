# pragma once

# include <string>

class Para
{
	public:
		Para(std::string filename = "");

		// file paths
		char fieldpath[1024];
		char probepath[1024];
		char statpath[1024];
		char postpath[1024];

		// physical parameters
		double Re;			// Reynolds number
		double init_ener;	// initial turbulent intensity ( = <u'^2 + v'^2 + w'^2> / 2 )
		
		// grid settings
		int Nx;				// grid number of streamwise direction
		int Ny;				// grid number of wall normal direction
		int Nz;				// grid number of spanwise direction
		double Lx;			// domain length of streamwise direction
		double Ly;			// channel height. Not adjustable (must be 2.0 if the non-dimensional equation is solved)
		double Lz;			// domain length of spanwise direction
		double dy_min;		// y coordinate at first grid off wall

		// time evolution control
		int Nt;				// total time steps to evolve
		double dt;			// time step length

		// output control
		int nwrite;			// step interval for writing whole fields (0 for not outputing)
		int nprint;			// step interval for printing and updating log of computation state (0 for not outputing)
		int nprobe;			// step interval for writing time serials of probed layers (0 for not outputing)

		// input control
		int nread;			// step number of whole field files to read for continuing computation (0 for not inputing)
		int *jprbs;			// layer indexes of layers to be probed (jprbs[0] store the number of layers)

		void readPara(std::string filename);	// if no filename if provided, initiate with debug parameters
		void showPara();

		int*	dim() { static int dim[3] = {Nx, Ny, Nz}; return (int*) dim; };
		double*	len() { static double len[3] = {Lx, Ly, Lz}; return (double*) len; };

	private:
		int* parseJprbs(char *str);
};






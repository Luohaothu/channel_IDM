# pragma once


class Matrix
{
	public:
		Matrix(int rank);
		~Matrix();

		void tdma (double *a, double *b, double *c, double *d);
		void ctdma(double *a, double *b, double *c, double *d);

	private:
		int n;
		// double *l;	// array l is not needed since it can be replaced by a single number l in execution
		double *u;
		double *rho;
		double *sig;

};



// NOTE:
// matrix functions are very important to the overall efficiency
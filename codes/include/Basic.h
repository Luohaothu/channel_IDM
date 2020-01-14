# pragma once

# include <fftw3.h>

# define PI 3.1415926535898




class Mesh
{
	public:
		Mesh(int n1, int n2, int n3, double l1, double l2, double l3, double dy_min=0, char* path=NULL);	// here n2=Ny is less than the array dimension by 1
		void freeall();
		void initYmesh();

		const int Nx, Ny, Nz, Nxz;
		const double Lx, Ly, Lz;		// domain lengths
		const double dx, dz, dx2, dz2;	// grid interval of X & Z and their squares
		const double vol;				// volume of the whole physical domain

		double *y, *yc;				// y coordinates in wall-normal direction (y for V, yc U,W,P)
		double *dy, *h;				// y interval of the wall-normal grid (dy for V, h for U,W,P)
		double *dvol;				// volume of every cell
		double *hm,  *hc,  *hp;		// coefficients for wall-normal 2nd-order derivative of U,W
		double *dym, *dyc, *dyp;	// coefficients for wall-normal 2nd-order derivative of V

		double *pmj, *pcj, *ppj;	// coefficients for wall-normal 2nd-order derivative of P in Poisson equation, boundaries excluded
		double *ak1, *ak3;			// coefficients of the Fourier transformed Laplacian operator (related to wavenumbers)

		int *kpa, *kma;				// the forward and backward node indices in periodic Z direction
		int *ipa, *ima;				// the forward and backward node indices in periodic X direction

		int IDX(int i, int j, int k) const { return Nxz * j + Nx * k + i; };
		double kx(int i) const { return ( i - ( i > (int)(Nx/2) ? Nx : 0 ) ) * (2.*PI/Lx); };
		double kz(int k) const { return ( k - ( k > (int)(Nz/2) ? Nz : 0 ) ) * (2.*PI/Lz); };

		void writeYmesh(char *path) const;

	private:
		void getYmesh(double dy_min); // generate hypertan y and infer yc
		void getYmesh(char *path);    // read existing y and infer yc

	// protected:
	// 	int idx, ip, im, jp, jm, kp, km,
	// 		ipjp, jpkp, ipkp,  ipjm, jpkm, imkp,
	// 		imjp, jmkp, ipkm,  imjm, jmkm, imkm;

	// 	void getIDX(int i, int j, int k)
	// 	{
	// 		idx = IDX(i,j,k);
	// 		ip = IDX(ipa[i],j,k); jp = IDX(i,j+1,k); kp = IDX(i,j,kpa[k]);
	// 		im = IDX(ima[i],j,k); jm = IDX(i,j-1,k); km = IDX(i,j,kma[k]);
	// 		ipjp = IDX(ipa[i],j+1,k); jpkp = IDX(i,j+1,kpa[k]); ipkp = IDX(ipa[i],j,kpa[k]);
	// 		ipjm = IDX(ipa[i],j-1,k); jpkm = IDX(i,j+1,kma[k]); imkp = IDX(ima[i],j,kpa[k]);
	// 		imjp = IDX(ima[i],j+1,k); jmkp = IDX(i,j-1,kpa[k]); ipkm = IDX(ipa[i],j,kma[k]);
	// 		imjm = IDX(ima[i],j-1,k); jmkm = IDX(i,j-1,kma[k]); imkm = IDX(ima[i],j,kma[k]);
	// 	};
};



class Bulk
{
	public:
		Bulk(int n1, int n2, int n3, bool inift=false);	// memory allocation, dim is the array dimension, unrelated to the mesh dimension
		~Bulk();

		double* fft ();
		double* ifft();

		// memory access functions
		double& id (int i, int j=0, int k=0) const { return  q[nxz * j + nx * k + i]; };
		double& idf(int i, int j=0, int k=0) const { return fq[nxzr* j + nxr* k + i]; };
		double* lyrGet (int j=0) const { return &( q[nxz * j]); };
		double* lyrGetF(int j=0) const { return &(fq[nxzr* j]); };
		double* blkGet () const { return  q; };
		double* blkGetF() const { return fq; };

		// tool functions
		Bulk& lyrSet(double a, int j=0) { this->layerTraverse(a, j, set); return *this; }; // set to a
		Bulk& lyrAdd(double a, int j=0) { this->layerTraverse(a, j, add); return *this; }; // add a
		Bulk& lyrMlt(double a, int j=0) { this->layerTraverse(a, j, mlt); return *this; }; // multiply by a
		Bulk& lyrSet(double *src, int j=0) { this->layerTraverse(src, j, set); return *this; }; // copy from src
		Bulk& lyrAdd(double *src, int j=0) { this->layerTraverse(src, j, add); return *this; }; // add src
		Bulk& lyrMns(double *src, int j=0) { this->layerTraverse(src, j, mns); return *this; }; // minus src
		Bulk& lyrMlt(double *src, int j=0) { this->layerTraverse(src, j, mlt); return *this; }; // multiply by src
		Bulk& lyrDvd(double *src, int j=0) { this->layerTraverse(src, j, dvd); return *this; }; // divide by src

		Bulk& blkSet(double a) { this->bulkTraverse(a, set); return *this; }; // note: bulk functions will change the boundary
		Bulk& blkAdd(double a) { this->bulkTraverse(a, add); return *this; };
		Bulk& blkMlt(double a) { this->bulkTraverse(a, mlt); return *this; };
		Bulk& blkSet(double *src) { this->bulkTraverse(src, set); return *this; };
		Bulk& blkAdd(double *src) { this->bulkTraverse(src, add); return *this; };
		Bulk& blkMns(double *src) { this->bulkTraverse(src, mns); return *this; };
		Bulk& blkMlt(double *src) { this->bulkTraverse(src, mlt); return *this; };
		Bulk& blkDvd(double *src) { this->bulkTraverse(src, dvd); return *this; };

		// IO functions
		void fileIO(char *path, char *name, char mode) const;
		void debug_AsciiOutput(char *path, char *name, int j1, int j2) const;

	protected:
		double *q;	// pointer to the bulk memory
		double *fq;	// Fourier tansformed q with only the kx >= 0 half

	private:
		const int nx, ny, nz, nxz;
		const int nxc, nxr, nxzc, nxzr;
		fftw_plan *frcs;	// list of fft plans: 2D FFT from real data to complex data, FORWARD (exponent -1) and no normalization
		fftw_plan *fcrs;	// list of ifft plans: 2D FFT from complex data to real data, BACKWARD (exponent 1) with normalization

		static void set(double &b, double a) {b = a; };
		static void add(double &b, double a) {b += a;};
		static void mns(double &b, double a) {b -= a;};
		static void mlt(double &b, double a) {b *= a;};
		static void dvd(double &b, double a) {b /= a;};

		void layerTraverse(double a,    int j, void (*pfun)(double &b, double a));
		void layerTraverse(double *src, int j, void (*pfun)(double &b, double a));
		void bulkTraverse (double a,           void (*pfun)(double &b, double a));
		void bulkTraverse (double *src,        void (*pfun)(double &b, double a));
};



class Scla: private Mesh, public Bulk
{
	public:
		Scla(const Mesh &mesh, bool inift=false):
			Mesh(mesh), Bulk(Nx, Ny+1, Nz, inift) {};

		Mesh& meshGet() const { return (Mesh&)(*this); };

		// integration
		double layerMean(int j=0) const;
		double yMeanU(int i, int k) const;
		double yMeanV(int i, int k) const;
		double bulkMeanU() const;
		double bulkMeanV() const;

		// interpolation
		void layerUG2CC(double *dst, int j);
		void layerVG2CC(double *dst, int j);
		void layerWG2CC(double *dst, int j);

		// operators
		double* gradient(int i, int j, int k);

		// double operator | (int j) const { return this->layerMean(j); };
		double av(int j=0) const { return this->layerMean(j); };
		double* operator [] (int j) const { return this->lyrGet(j); };

		Scla& operator = (double a) { this->blkSet(a); return *this; }; // note: bulk functions will change the boundary
		Scla& operator +=(double a) { this->blkAdd(a); return *this; };
		Scla& operator *=(double a) { this->blkMlt(a); return *this; };

		Scla& operator = (double *src) { this->blkSet(src); return *this; };
		Scla& operator +=(double *src) { this->blkAdd(src); return *this; };
		Scla& operator -=(double *src) { this->blkMns(src); return *this; };
		Scla& operator *=(double *src) { this->blkMlt(src); return *this; };
		Scla& operator /=(double *src) { this->blkDvd(src); return *this; };

		Scla& operator = (Scla &src) { this->blkSet(src.blkGet()); return *this; };
		Scla& operator +=(Scla &src) { this->blkAdd(src.blkGet()); return *this; };
		Scla& operator -=(Scla &src) { this->blkMns(src.blkGet()); return *this; };
		Scla& operator *=(Scla &src) { this->blkMlt(src.blkGet()); return *this; };
		Scla& operator /=(Scla &src) { this->blkDvd(src.blkGet()); return *this; };

		// Scla(const Mesh &mesh, const Bulk &bulk): Mesh(mesh), Bulk(bulk) {};
};


class Vctr: private Mesh
{
	public:
		Vctr(const Mesh &mesh): Mesh(mesh), com1(mesh), com2(mesh), com3(mesh)
		{ u = com1.blkGet(); v = com2.blkGet(); w = com3.blkGet(); };

		Scla& operator [] (int i) const { return (Scla&)(i==1 ? com1 : i==2 ? com2 : com3); };
		// Vctr& operator = (Vctr &src) { com1 = src[1]; com2 = src[2]; com3 = src[3]; return *this; };
		// Vctr& operator = (double a) { com1 = a; com2 = a; com3 = a; return *this; };

		Mesh& meshGet() const { return (Mesh&)(*this); };

		// vector interpolation
		void layer2CC(double *dst1, double *dst2, double *dst3, int j)
		{
			com1.layerUG2CC(dst1, j);
			com2.layerVG2CC(dst2, j);
			com3.layerWG2CC(dst3, j);
		};

		// vector operators
		double  divergence(int i, int j, int k) const;
		double  convection(int i, int j, int k) const;
		double* strainrate(int i, int j, int k) const;
		double* gradient  (int i, int j, int k) const;

	private:
		Scla com1, com2, com3;
		double *u, *v, *w;
		
		// Vctr(const Mesh &mesh, const Scla &s1, const Scla &s2, const Scla &s3):
		// 	Mesh(mesh),
		// 	com1(s1), com2(s2), com3(s3),
		// 	u(com1.bulkGet()), v(com2.bulkGet()), w(com3.bulkGet()) {};
};




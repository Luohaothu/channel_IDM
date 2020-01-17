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
};



class Bulk
{
	public:
		Bulk(int n1, int n2, int n3);	// memory allocation, dim is the array dimension, unrelated to the mesh dimension
		~Bulk();

		// memory access functions
		double& id (int idx) const { return q[idx]; };
		double& id (int i, int j, int k) const { return q[nxz * j + nx * k + i]; };
		double& idf(int i, int j, int k) const { return q[nxzr* j + nxr* k + i + nxz]; };
		double* lyrGet (int j=0) const { return &(q[nxz * j]); };
		double* lyrGetF(int j=0) const { return &(q[nxzr* j + nxz]); };
		double* blkGet () const { return lyrGet(); };
		double* blkGetF() const { return lyrGetF(); };

		// tool functions
		void fft ();
		void ifft();

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
		void fileIO(const char *path, const char *name, char mode) const;
		void debug_AsciiOutput(const char *path, const char *name, int j1, int j2, bool ifreal=true) const;

	private:
		double *q;	// pointer to the bulk memory

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
		Scla(const Mesh &mesh): Mesh(mesh), Bulk(Nx, Ny+1, Nz) {};

		const Mesh& meshGet() const { return (const Mesh&)(*this); };

		// integration
		double layerMean(int j=0) const;
		double yMeanU(int i, int k) const;
		double yMeanV(int i, int k) const;
		double bulkMeanU() const;
		double bulkMeanV() const;

		// interpolation
		void layerUG2CC(double *dst, int j) const;
		void layerVG2CC(double *dst, int j) const;
		void layerWG2CC(double *dst, int j) const;

		void layerCC2XE(double *dst, int j) const;
		void layerCC2YE(double *dst, int j) const;
		void layerCC2ZE(double *dst, int j) const;
		void CC2EG(double *dst1, double *dst2, double *dst3) const;

		// operators
		double* gradient(int i, int j, int k) const;

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
};


class Vctr: private Mesh
{
	public:
		Vctr(const Mesh &mesh): Mesh(mesh), V1(mesh), V2(mesh), V3(mesh) {};

		const Mesh& meshGet() const { return (const Mesh&)(*this); };

		Scla& operator [] (int i) const { return (Scla&)(i==1 ? V1 : i==2 ? V2 : V3); };

		void ptrGet(double *&u, double *&v, double *&w) const
		{ u = V1.blkGet(); v = V2.blkGet(); w = V3.blkGet(); };

		Vctr& reset(double a=0) { V1 = a; V2 = a; V3 = a; return (Vctr&)(*this); };

		// vector interpolation
		void layer2CC(double *dst1, double *dst2, double *dst3, int j) const
		{
			V1.layerUG2CC(dst1, j);
			V2.layerVG2CC(dst2, j);
			V3.layerWG2CC(dst3, j);
		};

		// vector operators
		double  divergence(int i, int j, int k) const;
		double  convection(int i, int j, int k) const;
		double* strainrate(int i, int j, int k) const;
		double* gradient  (int i, int j, int k) const;

	private:
		Scla V1, V2, V3;
};


class Feld: private Mesh
{
	public:
		Vctr V;
		Scla S;

		Feld(const Mesh &mesh): Mesh(mesh), V(mesh), S(mesh) {};

		const Mesh& meshGet() const { return (const Mesh&)(*this); };

		void ptrGet(double *&u, double *&v, double *&w, double *&p) const
		{ V.ptrGet(u, v, w); p = S.blkGet(); };

		Feld& reset(double a=0) { V.reset(a); S = a; return (Feld&)(*this); };
		Feld& initrand(double energy);
		Feld& initfrom(const Feld &feld);

		void CC2EG() { S.CC2EG(V[1].blkGet(), V[2].blkGet(), V[3].blkGet()); };

		// IO functions
		Feld& writeField(const char *path, int tstep, char *suffix) const;
		Feld& readField (const char *path, int tstep, char *suffix) const;
		void writeTecplot(const char *path, int tstep, double time) const;

	private:
		// Scla& operator [] (int i) const { return (Scla&)(i<4 ? V[i] : S) };
};



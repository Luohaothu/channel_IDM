#include "Basic.h"

using namespace std;


Mesh::Mesh(const Geometry &geo):
geo(geo),
Nx(geo.Nx),
Ny(geo.Ny),
Nz(geo.Nz),
Lx(geo.Lx),
Ly(geo.Ly),
Lz(geo.Lz),
Nxc(geo.Nxc),
Nzc(geo.Nzc),
Nxr(geo.Nxr),
Nzr(geo.Nzr)
{};


void Mesh::ipx(int i, int j, int k, int &ip, int &jp, int &kp) const
{
	ip = idx(ipa(i), j, k);
	jp = idx(i, jpa(j), k);
	kp = idx(i, j, kpa(k));
}
void Mesh::imx(int i, int j, int k, int &im, int &jm, int &km) const
{
	im = idx(ima(i), j, k);
	jm = idx(i, jma(j), k);
	km = idx(i, j, kma(k));
}
void Mesh::ppx(int i, int j, int k, int &ipjp, int &jpkp, int &ipkp) const
{
	int ip = ipa(i);
	int jp = jpa(j);
	int kp = kpa(k);
	ipjp = idx(ip, jp, k);
	jpkp = idx(i, jp, kp);
	ipkp = idx(ip, j, kp);
}
void Mesh::pmx(int i, int j, int k, int &ipjm, int &jpkm, int &ipkm) const
{
	int ip = ipa(i), jp = jpa(j);
	int jm = jma(j), km = kma(k);
	ipjm = idx(ip, jm, k);
	jpkm = idx(i, jp, km);
	ipkm = idx(ip, j, km);
}
void Mesh::mpx(int i, int j, int k, int &imjp, int &jmkp, int &imkp) const
{
	int im = ima(i), jm = jma(j);
	int jp = jpa(j), kp = kpa(k);
	imjp = idx(im, jp, k);
	jmkp = idx(i, jm, kp);
	imkp = idx(im, j, kp);
}
void Mesh::mmx(int i, int j, int k, int &imjm, int &jmkm, int &imkm) const
{
	int im = ima(i);
	int jm = jma(j);
	int km = kma(k);
	imjm = idx(im, jm, k);
	jmkm = idx(i, jm, km);
	imkm = idx(im, j, km);
}


void Mesh::dcx(int i, int j, int k, double &dxc, double &dyc, double &dzc) const
{
	dxc = dx(i);
	dyc = dy(j);
	dzc = dz(k);
}
void Mesh::dpx(int i, int j, int k, double &dxp, double &dyp, double &dzp) const
{
	dxp = dx(i+1);
	dyp = dy(j+1);
	dzp = dz(k+1);
}
void Mesh::dmx(int i, int j, int k, double &dxm, double &dym, double &dzm) const
{
	dxm = dx(i-1);
	dym = dy(j-1);
	dzm = dz(k-1);
}
void Mesh::hcx(int i, int j, int k, double &hxc, double &hyc, double &hzc) const
{
	hxc = hx(i);
	hyc = hy(j);
	hzc = hz(k);
}
void Mesh::hpx(int i, int j, int k, double &hxp, double &hyp, double &hzp) const
{
	hxp = hx(i+1);
	hyp = hy(j+1);
	hzp = hz(k+1);
}
void Mesh::hmx(int i, int j, int k, double &hxm, double &hym, double &hzm) const
{
	hxm = hx(i-1);
	hym = hy(j-1);
	hzm = hz(k-1);
}








// often used statements (for copy & paste)

// int id =              ms.idx(i,j,k);
// int ip, jp, kp;       ms.ipx(i,j,k,ip,jp,kp);
// int im, jm, km;       ms.imx(i,j,k,im,jm,km);
// int ipjp, jpkp, ipkp; ms.ppx(i,j,k,ipjp,jpkp,ipkp);
// int ipjm, jpkm, ipkm; ms.pmx(i,j,k,ipjm,jpkm,ipkm);
// int imjp, jmkp, imkp; ms.mpx(i,j,k,imjp,jmkp,imkp);
// int imjm, jmkm, imkm; ms.mmx(i,j,k,imjm,jmkm,imkm);

// double dxc, dyc, dzc; ms.dcx(i,j,k,dxc,dyc,dzc);
// double dxp, dyp, dzp; ms.dpx(i,j,k,dxp,dyp,dzp);
// double dxm, dym, dzm; ms.dmx(i,j,k,dxm,dym,dzm);
// double hxc, hyc, hzc; ms.hcx(i,j,k,hxc,hyc,hzc);
// double hxp, hyp, hzp; ms.hpx(i,j,k,hxp,hyp,hzp);
// double hxm, hym, hzm; ms.hmx(i,j,k,hxm,hym,hzm);

// for (int j=1; j<ms.Ny; j++) { double dym,dyp,dyc=ms.dy(j,dym,dyp), hym,hyp,hyc=ms.hy(j,hym,hyp);
// for (int k=1; k<ms.Nz; k++) { double dzm,dzp,dzc=ms.dz(k,dzm,dzp), hzm,hzp,hzc=ms.hz(k,hzm,hzp);
// for (int i=1; i<ms.Nx; i++) { double dxm,dxp,dxc=ms.dx(i,dxm,dxp), hxm,hxp,hxc=ms.hx(i,hxm,hxp);




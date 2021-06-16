import numpy as np
import os

import fileIO


class Write:
	def __init__(self, path, casename=None, irange=None, jrange=None, krange=None):
		self.workpath = path
		self.casename = casename
		self.irange = irange
		self.jrange = jrange
		self.krange = krange

	def wallscale(self, para):

		with open(self.workpath+"wallscale.txt", 'w') as fp:
			fp.write("Re_tau = %.18e\n"%para.Ret)
			fp.write("u_tau = %.18e\n"%para.utau)
			fp.write("tau_w = %.18e\n"%para.tauw)
			fp.write("delta_nu = %.18e\n"%para.dnu)
			fp.write("t_nu = %.18e\n"%para.tnu)
			fp.write("dy_min_plus = %.18e\n"%((para.y[2]-para.y[1])/para.dnu))
			fp.write("dy_max_plus = %.18e\n"%((para.y[para.Ny//2+1]-para.y[para.Ny//2])/para.dnu))

	def profiles(self, para, stas):

		header = \
			'Title = "profiles of basic statistics"\n' + \
			'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
			% (
				"y<sup>+</sup>",
				"<u><sup>+</sup>",
				"<v><sup>+</sup>",
				"<w><sup>+</sup>",
				"<p><sup>+</sup>",
				"<u'u'><sup>+</sup>",
				"<v'v'><sup>+</sup>",
				"<w'w'><sup>+</sup>",
				"<u'v'><sup>+</sup>",
				"<v'w'><sup>+</sup>",
				"<u'w'><sup>+</sup>",
				"<u'p'><sup>+</sup>",
				"<v'p'><sup>+</sup>",
				"<w'p'><sup>+</sup>",
				"<p'p'><sup>+</sup>",
				) + \
			'zone t = "%s", i = %i' %(self.casename, len(self.jrange))

		data = np.vstack([
			para.yc/para.lc,
			stas.Um/para.uc,
			stas.Vm/para.uc,
			stas.Wm/para.uc,
			stas.Pm/para.pc,
			stas.R11/para.uc**2,
			stas.R22/para.uc**2,
			stas.R33/para.uc**2,
			stas.R12/para.uc**2,
			stas.R23/para.uc**2,
			stas.R13/para.uc**2,
			stas.Rpu/para.uc/para.pc,
			stas.Rpv/para.uc/para.pc,
			stas.Rpw/para.uc/para.pc,
			stas.Rpp/para.pc**2,
			]).T[self.jrange]

		np.savetxt(self.workpath+"profiles.dat", data, header=header, comments='')

	def budgets(self, para, bgts):

		header = \
			'Title = "profiles of budgets"\n' + \
			'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
			% (
				"y<sup>+</sup>",
				"P<sup>+</sup>",
				"D<sub>t</sub><sup>+</sup>",
				"D<sub>p</sub><sup>+</sup>",
				"D<sub>v</sub><sup>+</sup>",
				"<greek>F</greek><sup>+</sup>",
				"<greek>e</greek><sup>+</sup>",
				"balance",
				) + \
			'zone t = "%s", i = %i' %(self.casename, len(self.jrange))

		data = np.vstack([
			para.yc/para.lc,
			bgts.Prod/(para.uc**3/para.lc),
			bgts.Dtur/(para.uc**3/para.lc),
			bgts.Dpre/(para.uc**3/para.lc),
			bgts.Dvis/(para.uc**3/para.lc),
			bgts.Rdis/(para.uc**3/para.lc),
			bgts.Epsl/(para.uc**3/para.lc),
			bgts.Baln/(para.uc**3/para.lc),
			]).T[self.jrange]

		np.savetxt(self.workpath+"budgets.dat", data, header=header, comments='')

	def es2d(self, para, stas):

		header = \
			'Title = "2D energy spectra"\n' + \
			'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
			% (
				"log<sub>10</sub>(<greek>l</greek><sub>x</sub><sup>+</sup>)",
				"log<sub>10</sub>(<greek>l</greek><sub>z</sub><sup>+</sup>)",
				"y<sup>+</sup>",
				"k<sub>x</sub>k<sub>z</sub>E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>k<sub>z</sub>E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>k<sub>z</sub>E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>k<sub>z</sub>E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>k<sub>z</sub>E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>k<sub>z</sub>E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>k<sub>z</sub>E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"
				) + \
			'zone t = "%s", i = %i, j = %i, k = %i' %(self.casename, len(self.irange), len(self.krange), len(self.jrange))

		data = np.array([
			np.empty_like(stas.Euu),
			np.empty_like(stas.Euu),
			np.empty_like(stas.Euu),
			stas.Euu / para.uc**2,
			stas.Evv / para.uc**2,
			stas.Eww / para.uc**2,
			stas.Epp / para.pc**2,
			stas.Euv / para.uc**2,
			stas.Evw / para.uc**2,
			stas.Euw / para.uc**2,
			])[:,self.jrange][:,:,self.krange][:,:,:,self.irange]

		irange_ = np.reshape(self.irange, [1,1,-1])
		krange_ = np.reshape(self.krange, [1,-1,1])
		jrange_ = np.reshape(self.jrange, [-1,1,1])

		# log wavenumber & inner-scaled
		data[0] = np.log10(2*np.pi/para.lc / para.kx[irange_])
		data[1] = np.log10(2*np.pi/para.lc / para.kz[krange_])
		data[2] = para.yc[jrange_]/para.lc

		# pre-multiply
		data[3:] *= para.kx[irange_] * para.kz[krange_] / (4*np.pi**2 / (para.Lx*para.Lz))

		np.savetxt(self.workpath+'es2d.dat', data.reshape([len(data),-1]).T, header=header, comments='')
		
		if not os.system("preplot " + self.workpath+'es2d.dat'):
			# os.system("rm -f " + pame)
			pass

	def es1d_xy(self, para, stas):

		header = \
			'Title = "1D streamwise energy spectra"\n' + \
			'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
			% (
				"log<sub>10</sub>(<greek>l</greek><sub>x</sub><sup>+</sup>)",
				"log<sub>10</sub>(y<sup>+</sup>)",
				"k<sub>x</sub>E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>x</sub>E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"
				) + \
			'zone t = "%s", i = %i, j = %i' %(self.casename, len(self.irange), len(self.jrange))


		data = np.array([
			np.empty_like(stas.Euu[:,0]),
			np.empty_like(stas.Euu[:,0]),
			np.sum(stas.Euu, axis=1) / para.uc**2,
			np.sum(stas.Evv, axis=1) / para.uc**2,
			np.sum(stas.Eww, axis=1) / para.uc**2,
			np.sum(stas.Epp, axis=1) / para.pc**2,
			np.sum(stas.Euv, axis=1) / para.uc**2,
			np.sum(stas.Evw, axis=1) / para.uc**2,
			np.sum(stas.Euw, axis=1) / para.uc**2,
			])[:,self.jrange][:,:,self.irange]

		irange_ = np.reshape(self.irange, [1,-1])
		jrange_ = np.reshape(self.jrange, [-1,1])

		# log wavenumber & inner-scaled
		data[0] = np.log10(2*np.pi/para.lc / para.kx[irange_])
		data[1] = np.log10(para.yc[jrange_]/para.lc)

		# pre-multiply
		data[2:] *= para.kx[irange_] / (2*np.pi/para.Lx)

		np.savetxt(self.workpath+'es1d_xy.dat', data.reshape([len(data),-1]).T, header=header, comments='')

		if not os.system("preplot " + self.workpath+'es1d_xy.dat'):
			# os.system("rm -f " + pame)
			pass

	def es1d_zy(self, para, stas):

		header = \
			'Title = "1D spanwise energy spectra"\n' + \
			'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
			% (
				"log<sub>10</sub>(<greek>l</greek><sub>z</sub><sup>+</sup>)",
				"log<sub>10</sub>(y<sup>+</sup>)",
				"k<sub>z</sub>E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>z</sub>E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>z</sub>E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>z</sub>E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>z</sub>E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>z</sub>E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
				"k<sub>z</sub>E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"
				) + \
			'zone t = "%s", i = %i, j = %i' %(self.casename, len(self.krange), len(self.jrange))

		data = np.array([
			np.empty_like(stas.Euu[:,0]),
			np.empty_like(stas.Euu[:,0]),
			np.sum(stas.Euu, axis=1) / para.uc**2,
			np.sum(stas.Evv, axis=1) / para.uc**2,
			np.sum(stas.Eww, axis=1) / para.uc**2,
			np.sum(stas.Epp, axis=1) / para.pc**2,
			np.sum(stas.Euv, axis=1) / para.uc**2,
			np.sum(stas.Evw, axis=1) / para.uc**2,
			np.sum(stas.Euw, axis=1) / para.uc**2,
			])[:,self.jrange][:,:,self.krange]

		krange_ = np.reshape(self.krange, [1,-1])
		jrange_ = np.reshape(self.jrange, [-1,1])

		# log wavenumber & inner-scaled
		data[0] = np.log10(2*np.pi/para.lc / para.kz[krange_])
		data[1] = np.log10(para.yc[jrange_]/para.lc)

		# pre-multiply
		data[2:] *= para.kz[krange_] / (2*np.pi/para.Lz)

		np.savetxt(self.workpath+'es1d_zy.dat', data.reshape([len(data),-1]).T, header=header, comments='')

		if not os.system("preplot " + self.workpath+'es1d_zy.dat'):
			# os.system("rm -f " + pame)
			pass

	def raw(self, para, stas):

		if 'raw' not in os.listdir(self.workpath): os.system('mkdir %sraw'%self.workpath)

		para.kx[:para.Nxc].astype(np.float64).tofile(self.workpath + 'raw/kxs.bin')
		para.kz.           astype(np.float64).tofile(self.workpath + 'raw/kzs.bin')
		para.yc[self.jrange].   astype(np.float64).tofile(self.workpath + 'raw/ys.bin')

		fileIO.write_channel(self.workpath + 'raw/Euu.bin', stas.Euu[self.jrange])
		fileIO.write_channel(self.workpath + 'raw/Evv.bin', stas.Evv[self.jrange])
		fileIO.write_channel(self.workpath + 'raw/Eww.bin', stas.Eww[self.jrange])
		fileIO.write_channel(self.workpath + 'raw/Epp.bin', stas.Epp[self.jrange])
		fileIO.write_channel(self.workpath + 'raw/Euvr.bin', stas.Euv.real[self.jrange])
		fileIO.write_channel(self.workpath + 'raw/Euvi.bin', stas.Euv.imag[self.jrange])
		fileIO.write_channel(self.workpath + 'raw/Evwr.bin', stas.Evw.real[self.jrange])
		fileIO.write_channel(self.workpath + 'raw/Evwi.bin', stas.Evw.imag[self.jrange])
		fileIO.write_channel(self.workpath + 'raw/Euwr.bin', stas.Euw.real[self.jrange])
		fileIO.write_channel(self.workpath + 'raw/Euwi.bin', stas.Euw.imag[self.jrange])




	def profiles_x(self, para, stas):

		header = \
			'Title = "profiles of basic statistics"\n' + \
			'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
			% (
				"x<sup>+</sup>",
				"y<sup>+</sup>",
				"U<sup>+</sup>",
				"V<sup>+</sup>",
				"W<sup>+</sup>",
				"P<sup>+</sup>",
				"<u'u'><sup>+</sup>",
				"<v'v'><sup>+</sup>",
				"<w'w'><sup>+</sup>",
				"<u'v'><sup>+</sup>",
				"<v'w'><sup>+</sup>",
				"<u'w'><sup>+</sup>",
				"<p'u'><sup>+</sup>",
				"<p'v'><sup>+</sup>",
				"<p'w'><sup>+</sup>",
				"<p'p'><sup>+</sup>"
				) + \
			'zone t = "%s", i = %i, j = %i' %(self.casename, len(self.irange), len(self.jrange))

		data = np.array([
			np.empty_like(stas.Um),
			np.empty_like(stas.Um),
			stas.Um /para.uc,
			stas.Vm /para.uc,
			stas.Wm /para.uc,
			stas.Pm /para.pc,
			stas.Ruu/para.uc**2,
			stas.Rvv/para.uc**2,
			stas.Rww/para.uc**2,
			stas.Ruv/para.uc**2,
			stas.Rvw/para.uc**2,
			stas.Ruw/para.uc**2,
			stas.Rpu/para.uc/para.pc,
			stas.Rpv/para.uc/para.pc,
			stas.Rpw/para.uc/para.pc,
			stas.Rpp/para.pc**2,
			])[:,self.jrange]

		data[:2] = np.meshgrid(
			para.xc[self.irange]/para.lc,
			para.yc[self.jrange]/para.lc )

		np.savetxt(self.workpath + 'profiles.dat', data.reshape([len(data), -1]).T, header=header, comments='')

		if not os.system("preplot " + self.workpath + 'profiles.dat'):
			# os.system("rm -f " + pame)
			pass

	def develops(self, para, stas):

		header = \
			'Title = "streamwise profiles of developing statistics"\n' + \
			'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
			% (
				"x<sup>+</sup>",
				"Re<sub><greek>q</greek></sub>",
				"Re<sub><greek>t</greek></sub>",
				"C<sub>f</sub>",
				"H<sub>12</sub>",
				"<greek>d</greek>",
				"<greek>d</greek><sup>*</sup>/<greek>d</greek>",
				"<greek>q</greek>/<greek>d</greek>",
				"<greek>d</greek><sup>**</sup>/<greek>d</greek>",
				) + \
			'zone t = "%s", i = %i' %(self.casename, len(self.irange))

		data = np.vstack([
			para.xc[1:-1]/para.lc,
			stas.Re_the,
			stas.Re_tau,
			stas.Cf,
			stas.H,
			stas.dlt /para.lc,
			stas.dlt1/para.lc,
			stas.dlt2/para.lc,
			stas.dlt3/para.lc,
			]).T[self.irange]

		np.savetxt(self.workpath + 'develops.dat', data, header=header, comments='')





import numpy as np
import os

from basic import DataSetInfo
import fileIO
from statis import Statis
from budgets import Budgets


para = DataSetInfo("mytest/")
para.scale_inner()

stas = Statis(para)
stas.calc_statis()
# stas.flipy()


bgts = Budgets(para)
bgts.calc_budgets()


with open(para.postpath+"wallscale.txt", 'w') as fp:
	fp.write("Re_tau = %.18e\n"%para.Ret)
	fp.write("u_tau = %.18e\n"%para.utau)
	fp.write("tau_w = %.18e\n"%para.tauw)
	fp.write("delta_nu = %.18e\n"%para.dnu)
	fp.write("t_nu = %.18e\n"%para.tnu)
	fp.write("dy_min_plus = %.18e\n"%((para.y[2]-para.y[1])/para.dnu))
	fp.write("dy_max_plus = %.18e\n"%((para.y[para.Ny//2+1]-para.y[para.Ny//2])/para.dnu))





casename = para.datapath.split('/')[-2]
jrange = range(1, para.Ny+1) # range(0, para.Ny//2+1) #
krange = range(1, para.Nzc)
irange = range(1, para.Nxc)




header = \
	'Title = "profiles of basic statistics"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"y<sup>+</sup>",
		"<u><sup>+</sup>", "<v><sup>+</sup>", "<w><sup>+</sup>", "<p><sup>+</sup>", \
		"<u'u'><sup>+</sup>", "<v'v'><sup>+</sup>", "<w'w'><sup>+</sup>", \
		"<u'v'><sup>+</sup>", "<v'w'><sup>+</sup>", "<u'w'><sup>+</sup>", \
		"<u'p'><sup>+</sup>", "<v'p'><sup>+</sup>", "<w'p'><sup>+</sup>", "<p'p'><sup>+</sup>"	) + \
	'zone t = "%s", i = %i' %( casename, len(jrange) )

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
	])

data = data.T[jrange]

np.savetxt(para.postpath+"profiles.dat", data, header=header, comments='')





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
	'zone t = "%s", i = %i' %( casename, len(jrange) )
data = np.vstack([
	para.yc/para.lc,
	bgts.Prod/(para.uc**3/para.lc),
	bgts.Dtur/(para.uc**3/para.lc),
	bgts.Dpre/(para.uc**3/para.lc),
	bgts.Dvis/(para.uc**3/para.lc),
	bgts.Rdis/(para.uc**3/para.lc),
	bgts.Epsl/(para.uc**3/para.lc),
	bgts.Baln/(para.uc**3/para.lc),
	])
data = data.T[jrange]
np.savetxt(para.postpath+"budgets.dat", data, header=header, comments='')





if 'raw' not in os.listdir(para.postpath): os.system('mkdir %sraw'%para.postpath)

para.kx[:para.Nxc].astype(np.float64).tofile(para.postpath + 'raw/kxs.bin')
para.kz.           astype(np.float64).tofile(para.postpath + 'raw/kzs.bin')
para.yc[jrange].   astype(np.float64).tofile(para.postpath + 'raw/ys.bin')

fileIO.write_channel(para.postpath + 'raw/Euu.bin', stas.Euu[jrange])
fileIO.write_channel(para.postpath + 'raw/Evv.bin', stas.Evv[jrange])
fileIO.write_channel(para.postpath + 'raw/Eww.bin', stas.Eww[jrange])
fileIO.write_channel(para.postpath + 'raw/Epp.bin', stas.Epp[jrange])
fileIO.write_channel(para.postpath + 'raw/Euvr.bin', stas.Euv.real[jrange])
fileIO.write_channel(para.postpath + 'raw/Euvi.bin', stas.Euv.imag[jrange])
fileIO.write_channel(para.postpath + 'raw/Evwr.bin', stas.Evw.real[jrange])
fileIO.write_channel(para.postpath + 'raw/Evwi.bin', stas.Evw.imag[jrange])
fileIO.write_channel(para.postpath + 'raw/Euwr.bin', stas.Euw.real[jrange])
fileIO.write_channel(para.postpath + 'raw/Euwi.bin', stas.Euw.imag[jrange])

stas.flipk()





header = \
	'Title = "2D energy spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"log<sub>10</sub>(<greek>l</greek><sub>x</sub><sup>+</sup>)",
		"log<sub>10</sub>(<greek>l</greek><sub>z</sub><sup>+</sup>)",
		"y<sup>+</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>k<sub>z</sub>E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"	) + \
	'zone t = "%s", i = %i, j = %i, k = %i' %( casename, len(jrange), len(krange), len(irange) )

data = np.empty([10, len(irange), len(krange), len(jrange)])
for i in irange:
	for k in krange:
		for j in jrange:
			data[:, irange.index(i), krange.index(k), jrange.index(j)] = [
				para.kx[i], para.kz[k], para.yc[j],
				stas.Euu[j,k,i],
				stas.Evv[j,k,i],
				stas.Eww[j,k,i],
				stas.Epp[j,k,i],
				stas.Euv[j,k,i],
				stas.Evw[j,k,i],
				stas.Euw[j,k,i]	]

data[3:] *= data[0] * data[1] / (4*np.pi**2 / para.Lx / para.Lz) / para.uc**2
data[6] *= para.uc**2 / para.pc**2
data[:2] = np.log10(2*np.pi / data[:2] / para.lc)
data[2] /= para.lc
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES2D.dat"
np.savetxt(pame, data, header=header, comments='')
if not os.system("preplot " + pame):
	# os.system("rm -f " + pame)
	pass





header = \
	'Title = "1D streamwise energy spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"log<sub>10</sub>(<greek>l</greek><sub>x</sub><sup>+</sup>)",
		"log<sub>10</sub>(y<sup>+</sup>)",
		"k<sub>x</sub>E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>x</sub>E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"	) + \
	'zone t = "%s", i = %i, j = %i' %( casename, len(jrange), len(irange) )

data = np.empty([9, len(irange), len(jrange)])
for i in irange:
	for j in jrange:
		data[:, irange.index(i), jrange.index(j)] = [
			para.kx[i], para.yc[j],
			np.sum(stas.Euu[j,:,i]),
			np.sum(stas.Evv[j,:,i]),
			np.sum(stas.Eww[j,:,i]),
			np.sum(stas.Epp[j,:,i]),
			np.sum(stas.Euv[j,:,i]),
			np.sum(stas.Evw[j,:,i]),
			np.sum(stas.Euw[j,:,i])	]

data[2:] *= data[0] / (2*np.pi / para.Lx) / para.uc**2
data[5] *= para.uc**2 / para.pc**2
data[0] = 2*np.pi / data[0]
data[:2] = np.log10(data[:2] / para.lc)
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES1D_xy.dat"
np.savetxt(pame, data, header=header, comments='')
if not os.system("preplot " + pame):
	# os.system("rm -f " + pame)
	pass





header = \
	'Title = "1D spanwise energy spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"log<sub>10</sub>(<greek>l</greek><sub>z</sub><sup>+</sup>)",
		"log<sub>10</sub>(y<sup>+</sup>)",
		"k<sub>z</sub>E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"k<sub>z</sub>E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"	) + \
	'zone t = "%s", i = %i, j = %i' %( casename, len(jrange), len(krange) )

data = np.empty([9, len(krange), len(jrange)])
for k in krange:
	for j in jrange:
		data[:, krange.index(k), jrange.index(j)] = [
			para.kz[k], para.yc[j],
			np.sum(stas.Euu[j,k]),
			np.sum(stas.Evv[j,k]),
			np.sum(stas.Eww[j,k]),
			np.sum(stas.Epp[j,k]),
			np.sum(stas.Euv[j,k]),
			np.sum(stas.Evw[j,k]),
			np.sum(stas.Euw[j,k])	]

data[2:] *= data[0] / (2*np.pi / para.Lz) / para.uc**2
data[5] *= para.uc**2 / para.pc**2
data[0] = 2*np.pi / data[0]
data[:2] = np.log10(data[:2] / para.lc)
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES1D_zy.dat"
np.savetxt(pame, data, header=header, comments='')
if not os.system("preplot " + pame):
	# os.system("rm -f " + pame)
	pass






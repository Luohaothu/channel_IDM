#!/root/Software/anaconda3/bin/python3
from basic import *
from statis import Statis
from budgets import Budgets


para = DataSetInfo("../data_test/test/")

stas = Statis(Field(para))
# bgts = Budgets(para)

stas.calc_statis()
stas.calc_wallscale()
stas.flipy()

# bgts.dissipation()
# bgts.flipy()



with open(para.postpath+"wallscale.txt", 'w') as fp:
	fp.write("Re_tau = %.18e\n"%stas.Ret)
	fp.write("u_tau = %.18e\n"%stas.utau)
	fp.write("tau_w = %.18e\n"%stas.tauw)
	fp.write("delta_nu = %.18e\n"%stas.dnu)
	fp.write("t_nu = %.18e\n"%stas.tnu)
	fp.write("dy_min_plus = %.18e\n"%((para.y[2]-para.y[1])/stas.dnu))
	fp.write("dy_max_plus = %.18e\n"%((para.y[para.Ny//2+1]-para.y[para.Ny//2])/stas.dnu))


qs = { 'Euu.bin':stas.Euu, 'Evv.bin':stas.Evv, 'Eww.bin':stas.Eww, 'Epp.bin':stas.Epp,
       'Euv.bin':stas.Euv, 'Evw.bin':stas.Evw, 'Euw.bin':stas.Euw, }
##################################################
for name in qs: write_channel(para.postpath + name, qs[name])
########### above: write & below: read ###########
# for name in qs: read_channel (para.postpath + name)
##################################################




# stas.outer_scale()
stas.inner_scale()

casename = para.datapath.split('/')[-2]
jrange = range(0, para.Ny//2+1) #range(1, para.Ny) #
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

data = np.vstack([ para.yc/stas.lc,
	stas.Um,	stas.Vm,	stas.Wm,	stas.Pm/stas.pc,
	stas.R11,	stas.R22,	stas.R33,	stas.R12,	stas.R23,	stas.R13,
	stas.Rpu,	stas.Rpv,	stas.Rpw,	stas.Rpp/stas.pc**2	])
data[1:4] /= stas.uc
data[5:11] /= stas.uc**2
data[11:14] /= stas.uc * stas.pc
data = data.T[jrange]

np.savetxt(para.postpath+"profiles.dat", data, header=header, comments='')





# header = \
# 	'Title = "profiles of budgets"\n' + \
# 	'variables = "%s", "%s"\n' % ( "y<sup>+</sup>", "<greek>e</greek><sup>+</sup>" ) + \
# 	'zone t = "%s", i = %i' %( casename, len(jrange) )
# data = np.vstack([ para.yc/stas.lc, bgts.epsl/(stas.uc**3/stas.lc) ])
# data = data.T[jrange]
# np.savetxt(para.postpath+"budgets.dat", data, header=header, comments='')





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

data[3:] *= data[0] * data[1] / (4*np.pi**2 / para.Lx / para.Lz) / stas.uc**2
data[6] *= stas.uc**2 / stas.pc**2
data[:2] = np.log10(2*np.pi / data[:2] / stas.lc)
data[2] /= stas.lc
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES2D.dat"
np.savetxt(pame, data, header=header, comments='')
if not system("preplot " + pame):
	system("rm -f " + pame)





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

data[2:] *= data[0] / (2*np.pi / para.Lx) / stas.uc**2
data[5] *= stas.uc**2 / stas.pc**2
data[0] = 2*np.pi / data[0]
data[:2] = np.log10(data[:2] / stas.lc)
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES1D_xy.dat"
np.savetxt(pame, data, header=header, comments='')
if not system("preplot " + pame):
	system("rm -f " + pame)





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

data[2:] *= data[0] / (2*np.pi / para.Lz) / stas.uc**2
data[5] *= stas.uc**2 / stas.pc**2
data[0] = 2*np.pi / data[0]
data[:2] = np.log10(data[:2] / stas.lc)
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES1D_zy.dat"
np.savetxt(pame, data, header=header, comments='')
if not system("preplot " + pame):
	system("rm -f " + pame)






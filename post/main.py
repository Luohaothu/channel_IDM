#!/root/Software/anaconda3/bin/python3
from basic import *

para = DataSetInfo("/run/media/student/DATA/whn/channel_IDM/data_LES1000_SM_CORS/")
feld = Field(para)
stas = Statis(para, feld)

stas.calc_statis()
stas.inner_scale()
stas.flipy()


qs = (stas.Euu, stas.Evv, stas.Eww, stas.Epp, stas.Euv, stas.Evw, stas.Euw)
fns= ("Euu.bin","Evv.bin","Eww.bin","Epp.bin","Euv.bin","Evw.bin","Euw.bin")
for q, fn in zip(qs, fns): write_channel(para.postpath+fn, q)

# nx, ny, nz = para.Nx, para.Ny + 1, para.Nz
# nxc, nzc = int(nx/2+1), int(nz/2+1)
# stas.Euu, stas.Evv, stas.Eww, stas.Epp= [np.zeros([ny, nzc, nxc]) for n in range(4)]
# stas.Euv, stas.Evw, stas.Euw = [np.zeros([ny, nzc, nxc]) for n in range(3)]

# qs = (stas.Euu, stas.Evv, stas.Eww, stas.Epp, stas.Euv, stas.Evw, stas.Euw)
# fns= ("Euu.bin","Evv.bin","Eww.bin","Epp.bin","Euv.bin","Evw.bin","Euw.bin")
# for q, fn in zip(qs, fns): q[:] = read_channel(para.postpath+fn)

with open(para.postpath+"innerscale.txt", 'w') as fp:
	fp.write("Re_tau = %.18e\n"%stas.Ret)
	fp.write("u_tau = %.18e\n"%stas.utau)
	fp.write("tau_w = %.18e\n"%stas.tauw)
	fp.write("delta_nu = %.18e\n"%stas.dnu)
	fp.write("t_nu = %.18e\n"%stas.tnu)
	fp.write("dy_min_plus = %.18e\n"%((para.y[2]-para.y[1])/stas.dnu))
	fp.write("dy_max_plus = %.18e\n"%((para.y[int(para.Ny/2+1)]-para.y[int(para.Ny/2)])/stas.dnu))




casename = para.datapath.split('/')[-2]
jrange = range(1, int(para.Ny/2+1)) #range(1, para.Ny) #
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

data = np.vstack([ para.yc/stas.dnu,
	stas.Um,	stas.Vm,	stas.Wm,	stas.Pm/stas.tauw,
	stas.R11,	stas.R22,	stas.R33,	stas.R12,	stas.R23,	stas.R13,
	stas.Rpu,	stas.Rpv,	stas.Rpw,	stas.Rpp/stas.tauw**2	])
data[1:4] /= stas.utau
data[5:11] /= stas.utau**2
data[11:14] /= stas.utau * stas.tauw
data = data.T[jrange]

np.savetxt(para.postpath+"profiles.dat", data, header=header, comments='')





header = \
	'Title = "2D energy spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"log<sub>10</sub>(<greek>l</greek><sub>x</sub><sup>+</sup>)",
		"log<sub>10</sub>(<greek>l</greek><sub>z</sub><sup>+</sup>)",
		"y<sup>+</sup>",
		"E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"	) + \
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

data[3:] *= data[0] * data[1] / (4*np.pi**2 / para.Lx / para.Lz) / stas.utau**2
data[6] *= stas.utau**2 / stas.tauw**2
data[:2] = np.log10(2*np.pi / data[:2] / stas.dnu)
data[2] /= stas.dnu
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES2D.dat"
np.savetxt(pame, data, header=header, comments='')
os.system("preplot " + pame)
os.system("rm -f " + pame)





header = \
	'Title = "1D streamwise energy spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"log<sub>10</sub>(<greek>l</greek><sub>x</sub><sup>+</sup>)",
		"log<sub>10</sub>(y<sup>+</sup>)",
		"E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"	) + \
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

data[2:] *= data[0] / (2*np.pi / para.Lx) / stas.utau**2
data[5] *= stas.utau**2 / stas.tauw**2
data[0] = 2*np.pi / data[0]
data[:2] = np.log10(data[:2] / stas.dnu)
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES1D_xy.dat"
np.savetxt(pame, data, header=header, comments='')
os.system("preplot " + pame)
os.system("rm -f " + pame)





header = \
	'Title = "1D spanwise energy spectra"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"log<sub>10</sub>(<greek>l</greek><sub>z</sub><sup>+</sup>)",
		"log<sub>10</sub>(y<sup>+</sup>)",
		"E<sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>vv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>ww</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>pp</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>uv</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>vw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>",
		"E<sub>uw</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"	) + \
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

data[2:] *= data[0] / (2*np.pi / para.Lz) / stas.utau**2
data[5] *= stas.utau**2 / stas.tauw**2
data[0] = 2*np.pi / data[0]
data[:2] = np.log10(data[:2] / stas.dnu)
data = np.array([np.ravel(temp) for temp in data]).T

pame = para.postpath + "ES1D_zy.dat"
np.savetxt(pame, data, header=header, comments='')
os.system("preplot " + pame)
os.system("rm -f " + pame)






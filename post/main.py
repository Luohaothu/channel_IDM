#!/root/Software/anaconda3/bin/python3
from basic import *

para = DataSetInfo("/run/media/student/DATA/whn/channel_IDM/SR_VEL/")
feld = Field(para)
stas = Statis(para, feld)

stas.calc_profs()
stas.inner_scale()

casename = para.datapath.split('/')[-2]
jrange = range(1, int(para.Ny/2)+1)


header = \
	'Title = "profiles of basic statistics"\n' + \
	'variables = "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"\n' \
	% (	"y<sup>+</sup>",
		"U<sup>+</sup>", "V<sup>+</sup>", "W<sup>+</sup>", "P<sup>+</sup>", \
		"<u'u'><sup>+</sup>", "<v'v'><sup>+</sup>", "<w'w'><sup>+</sup>", \
		"<u'v'><sup>+</sup>", "<v'w'><sup>+</sup>", "<u'w'><sup>+</sup>", \
		"<u'p'><sup>+</sup>", "<v'p'><sup>+</sup>", "<w'p'><sup>+</sup>", \
		"<p'p'><sup>+</sup>"	) + \
	'zone t = "%s", i = %i' \
	%(	casename, len(jrange)	)

data = np.vstack([ para.yc/stas.dnu,
	stas.Um,	stas.Vm,	stas.Wm,	stas.Pm/stas.tauw,
	stas.R11,	stas.R22,	stas.R33,	stas.R12,	stas.R23,	stas.R13,
	stas.Rpu,	stas.Rpv,	stas.Rpw,	stas.Rpp/stas.tauw**2	])
data[1:4] /= stas.utau
data[5:11] /= stas.utau**2
data[11:14] /= stas.utau * stas.tauw
data = data.T[jrange]

np.savetxt(para.postpath+"profiles.dat", data, header=header, comments='')






# r"$y^+$",
# r"$U^+$", r"$V^+$", r"$W^+$", r"$P^+$",
# r"$<u'u'>^+$", r"$<v'v'>^+$", r"$<w'w'>^+$",
# r"$<u'v'>^+$", r"$<v'w'>^+$", r"$<u'w'>^+$",
# r"$<u'p'>^+$", r"$<v'p'>^+$", r"$<w'p'>^+$", r"$<p'p'>^+$"
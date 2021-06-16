import numpy as np
import os

from basic import DataSetInfo
from statis import Statis
from budgets import Budgets
from write_statis import Write


if __name__ == '__main__':

	para = DataSetInfo("/mnt/disk2/whn/etbl/EBL_1420_small/")
	para.scale_inner()

	stas = Statis(para)
	stas.calc_statis()
	# stas.flipy()

	# bgts = Budgets(para)
	# bgts.calc_budgets()


	writer = Write(para.postpath,
		casename = para.datapath.split('/')[-2],
		irange = range(1, para.Nxc),
		krange = range(1, para.Nzc),
		jrange = range(1, para.Ny+1) # range(0, para.Ny//2+1) #,
		)

	writer.wallscale(para)
	writer.profiles(para, stas)
	writer.raw(para, stas)

	stas.flipk()
	writer.es2d(para, stas)
	writer.es1d_xy(para, stas)
	writer.es1d_zy(para, stas)

	exit()








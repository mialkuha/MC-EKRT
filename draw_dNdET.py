import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

#y_walls = np.linspace(-7.25, 7.25, 30)
et_walls = [ 5.45664, 7.27108, 9.68884, 12.9106, 17.2036, 22.9241, 30.5467, 40.7041, 54.2389, 72.2744, 96.3069, 128.331, 171.003, 227.865, 303.634, 404.598, 539.134, 718.406, 957.289, 1275.6, 1699.77, 2264.97, 3018.11, 4021.69, 5358.98, 7140.94, 9515.43, 12679.5, 16895.6, 22513.8 ]
#y_walls = [-7.5175, -7.13199, -6.74647, -6.36096, -5.97545, -5.58994, -5.20442, -4.81891, -4.4334, -4.04788, -3.66237, -3.27686, -2.89135, -2.50583, -2.12032, -1.73481, -1.34929, -0.963782, -0.578269, -0.192756, 0.192756, 0.578269, 0.963782, 1.34929, 1.73481, 2.12032, 2.50583, 2.89135, 3.27686, 3.66237, 4.04788, 4.4334, 4.81891, 5.20442, 5.58994, 5.97545, 6.36096, 6.74647, 7.13199, 7.5175]
#y_points = [y1+(y2-y1)/2 for (y1,y2) in zip(y_walls, y_walls[1:])]

obs_names = ['dNdET']

for name in obs_names:
	infile_names = [ 'dNdET_sim_testAA_mb.dat',
			 	     'dNdET_sim_testAA_mb_MC.dat',
			 	     'dNdET_sim_testAA_mb_MC_ND.dat',
			 	     'dNdET_sim_testAA_mb_SAT.dat',
			 	     'dNdET_sim_testAA_mb_SAT_MC.dat']

	pdf_name  = name+'.pdf'

	infiles = [ open(n,'r') for n in infile_names ]

	liness = [ f.readlines() for f in infiles ]

	for f in infiles:
		f.close()
		
	datas = [[float(c) for c in lines[4].split()] for lines in liness]

	# use latex for font rendering
	mpl.rcParams['text.usetex'] = True
	mpl.rcParams['legend.edgecolor'] = 'white'  

	fig, ax = plt.subplots()
	# plt.title(r'$\textrm{EPPS16, CT14Lo}$', size='xx-large')
	#ax.set_xlabel(r'$x$', size='xx-large')
	#ax.set_ylabel(r'$R_g * f_g$', size='xx-large')
	#plt.ylim(0,100)
	#ax.set_yscale('log')
	ax.set_xscale('log')
	plt.xlim(5.4,10000)
	plt.grid(visible=True, which='both', axis='both')
	ax.plot(et_walls, np.transpose(datas[0]), 'b-', label=r'$NOTHING$', linewidth=1)
	ax.plot(et_walls, np.transpose(datas[1]), 'r-', label=r'$MC$', linewidth=1)
	ax.plot(et_walls, np.transpose(datas[2]), 'g-', label=r'$MC\;ND$', linewidth=1)
	ax.plot(et_walls, np.transpose(datas[3]), 'b--', label=r'$SAT$', linewidth=1)
	ax.plot(et_walls, np.transpose(datas[4]), 'r--', label=r'$SAT\;MC$', linewidth=1)
	ax.minorticks_on()
	#ax.set_xticks((0.0001, 0.001, 0.01, 0.1, 1.0))
	#ax.set_xticklabels(('$1e-4$', '$1e-3$', '$1e-2$', '$0.1$', '$1.0$'), size='xx-large')
	#ax.set_yticks((0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8))
	#ax.set_yticklabels(('$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$', '$1.2$', '$1.4$', '$1.6$', '$1.8$'), size='xx-large')
	ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
	plt.legend(loc='upper right', ncol=1, fontsize='x-large')
	pp = PdfPages(pdf_name)
	plt.savefig(pp, format='pdf')
	pp.close()

	print('Printed to ' + pdf_name)

#b_walls = np.linspace(0, 20, 21)
#b2_walls = [b*b for b in b_walls]
#b_points = [b1+(b2-b1)/2 for (b1,b2) in zip(b_walls, b_walls[1:])]
#b2_points = [b1+(b2-b1)/2 for (b1,b2) in zip(b2_walls, b2_walls[1:])]
#
#obs_names = ['dETdb','dEdb']
#
#for name in obs_names:
#	infile_name1  = name + '_sim_pp_100k.dat'
#	infile_name2  = name + '_sim_AA_100k.dat'
#	infile_name3  = name + '_sim_sAA_100k.dat'
#
#	pdf_name  = name+'.pdf'
#
#	infile1  = open(infile_name1,'r')
#	infile2  = open(infile_name2,'r')
#	infile3  = open(infile_name3,'r')
#
#	lines1 = infile1.readlines()	
#	infile1.close()
#	lines2 = infile2.readlines()	
#	infile2.close()
#	lines3 = infile3.readlines()	
#	infile3.close()
#
#	data1  = []
#	data2  = []
#	data3  = []
#
#	for line in lines1:
#		if line.startswith('///') or len(line) < 4 :
#			continue
#		data1.append(list(map(float,line.split())))
#	for line in lines2:
#		if line.startswith('///') or len(line) < 4 :
#			continue
#		data2.append(list(map(float,line.split())))
#	for line in lines3:
#		if line.startswith('///') or len(line) < 4 :
#			continue
#		data3.append(list(map(float,line.split())))
#
#	# use latex for font rendering
#	mpl.rcParams['text.usetex'] = True
#	mpl.rcParams['legend.edgecolor'] = 'white'  
#
#	fig, ax = plt.subplots()
#	# plt.title(r'$\textrm{EPPS16, CT14Lo}$', size='xx-large')
#	#ax.set_xlabel(r'$x$', size='xx-large')
#	#ax.set_ylabel(r'$R_g * f_g$', size='xx-large')
#	#plt.ylim(0,10)
#	plt.xlim(0,20)
#	plt.grid(visible=True, which='both', axis='both')
#	ax.plot(y_walls, np.transpose(data1), 'b-', label=r'$pp$', linewidth=1)
#	ax.plot(y_walls, np.transpose(data2), 'r-', label=r'$AA$', linewidth=1)
#	ax.plot(y_walls, np.transpose(data3), 'g-', label=r'$sAA$', linewidth=1)
#	ax.minorticks_on()
#	#ax.set_xticks((0.0001, 0.001, 0.01, 0.1, 1.0))
#	#ax.set_xticklabels(('$1e-4$', '$1e-3$', '$1e-2$', '$0.1$', '$1.0$'), size='xx-large')
#	#ax.set_yticks((0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8))
#	#ax.set_yticklabels(('$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$', '$1.2$', '$1.4$', '$1.6$', '$1.8$'), size='xx-large')
#	ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
#	plt.legend(loc='upper right', ncol=1, fontsize='x-large')
#	pp = PdfPages(pdf_name)
#	plt.savefig(pp, format='pdf')
#	pp.close()
#
#	print('Printed to ' + pdf_name)
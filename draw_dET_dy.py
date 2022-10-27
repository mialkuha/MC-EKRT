import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

#y_walls = np.linspace(-7.25, 7.25, 30)
y_walls = [-7.5175, -7.13199, -6.74647, -6.36096, -5.97545, -5.58994, -5.20442, -4.81891, -4.4334, -4.04788, -3.66237, -3.27686, -2.89135, -2.50583, -2.12032, -1.73481, -1.34929, -0.963782, -0.578269, -0.192756, 0.192756, 0.578269, 0.963782, 1.34929, 1.73481, 2.12032, 2.50583, 2.89135, 3.27686, 3.66237, 4.04788, 4.4334, 4.81891, 5.20442, 5.58994, 5.97545, 6.36096, 6.74647, 7.13199, 7.5175]
y_points = [y1+(y2-y1)/2 for (y1,y2) in zip(y_walls, y_walls[1:])]

obs_names = ['dETdy','dEdy']

for name in obs_names:
	#infile_names = [ name + '_0_5_ekrt_like_variant_K=6.dat',
	#		 	     name + '_0_5_ekrt_like_variant_K=6_nPDF.dat',
	#		 	     name + '_0_5_ekrt_like_variant_K=6_SAT.dat',
	#		 	     name + '_0_5_ekrt_like_variant_K=6_nPDF_SAT.dat',
    #                 name + '_0_5_ekrt_like_variant_K=6_SAT_MC.dat',
    #                 name + '_0_5_ekrt_like_variant_K=6_nPDF_SAT_MC.dat']
	infile_names = [ name + '_0_5_ekrt_like_K=1_M=5.8_SAT.dat',
			 	     name + '_0_5_ekrt_like_K=2_M=2.4_SAT.dat',
			 	     name + '_0_5_ekrt_like_K=2_M=4_nPDF_SAT.dat',
                     name + '_0_5_ekrt_like_K=1_M=5.8_SAT_MC.dat',
			 	     name + '_0_5_ekrt_like_K=2_M=2.4_SAT_MC.dat',
			 	     name + '_0_5_ekrt_like_K=2_M=4_nPDF_SAT_MC.dat']

	pdf_name  = name+'_tau_order_comparison+MC_fixd.pdf'

	infiles = [ open(n,'r') for n in infile_names ]

	liness = [ f.readlines() for f in infiles ]

	for f in infiles:
		f.close()
		
	datas = [[float(c) for c in lines[4].split()] for lines in liness]

	# use latex for font rendering
	#mpl.rcParams['text.usetex'] = True
	mpl.rcParams['legend.edgecolor'] = 'white'  

	fig, ax = plt.subplots()
	# plt.title(r'$\textrm{EPPS16, CT14Lo}$', size='xx-large')
	#ax.set_xlabel(r'$x$', size='xx-large')
	#ax.set_ylabel(r'$R_g * f_g$', size='xx-large')
	#plt.ylim(0,10)
	plt.xlim(-8,8)
	plt.grid(visible=True, which='both', axis='both')
	#ax.plot(y_points, np.transpose(datas[0]), 'b-.', label='PDF', linewidth=1)
	#ax.plot(y_points, np.transpose(datas[1]), 'r-.', label='nPDF', linewidth=1)
	#ax.plot(y_points, np.transpose(datas[2]), 'b-', label='PDF SAT', linewidth=1)
	#ax.plot(y_points, np.transpose(datas[3]), 'r-', label='nPDF SAT', linewidth=1)
	#ax.plot(y_points, np.transpose(datas[4]), 'b--', label='PDF SAT+MC', linewidth=1)
	#ax.plot(y_points, np.transpose(datas[5]), 'r--', label='nPDF SAT+MC', linewidth=1)
	ax.plot(y_points, np.transpose(datas[0]), 'b-', label='K=1 M=5.8 PDF', linewidth=1)
	ax.plot(y_points, np.transpose(datas[1]), 'r-', label='K=2 M=2.4 PDF', linewidth=1)
	ax.plot(y_points, np.transpose(datas[2]), 'g-', label='K=2 M=4 nPDF', linewidth=1)
	ax.plot(y_points, np.transpose(datas[3]), 'b--', label='MC K=1 M=5.8 PDF', linewidth=1)
	ax.plot(y_points, np.transpose(datas[4]), 'r--', label='MC K=2 M=2.4 PDF', linewidth=1)
	ax.plot(y_points, np.transpose(datas[5]), 'g--', label='MC K=2 M=4 nPDF', linewidth=1)
	ax.minorticks_on()
	#ax.set_xticks((0.0001, 0.001, 0.01, 0.1, 1.0))
	#ax.set_xticklabels(('$1e-4$', '$1e-3$', '$1e-2$', '$0.1$', '$1.0$'), size='xx-large')
	#ax.set_yticks((0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8))
	#ax.set_yticklabels(('$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$', '$1.2$', '$1.4$', '$1.6$', '$1.8$'), size='xx-large')
	ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
	plt.legend(loc='lower right', ncol=1, fontsize='large')
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

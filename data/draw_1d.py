import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

#y_walls = np.linspace(-7.25, 7.25, 30)
y_walls = [-7.5175, -7.13199, -6.74647, -6.36096, -5.97545, -5.58994, -5.20442, -4.81891, -4.4334, -4.04788, -3.66237, -3.27686, -2.89135, -2.50583, -2.12032, -1.73481, -1.34929, -0.963782, -0.578269, -0.192756, 0.192756, 0.578269, 0.963782, 1.34929, 1.73481, 2.12032, 2.50583, 2.89135, 3.27686, 3.66237, 4.04788, 4.4334, 4.81891, 5.20442, 5.58994, 5.97545, 6.36096, 6.74647, 7.13199, 7.5175]
y_points = [y1+(y2-y1)/2 for (y1,y2) in zip(y_walls, y_walls[1:])]

obs_names = ['dETdy','dEdy','dNdy']

for name in obs_names:
	infile_name1  = name + '_sim_pp_100k.dat'
	infile_name2  = name + '_sim_AA_100k.dat'
	infile_name3  = name + '_sim_sAA_100k.dat'
	infile_name4  = name + '_sim_pp_100k_MC.dat'
	infile_name5  = name + '_sim_AA_100k_MC.dat'
	infile_name6  = name + '_sim_sAA_100k_MC.dat'

	pdf_name  = name+'.pdf'

	infile1  = open(infile_name1,'r')
	infile2  = open(infile_name2,'r')
	infile3  = open(infile_name3,'r')
	infile4  = open(infile_name4,'r')
	infile5  = open(infile_name5,'r')
	infile6  = open(infile_name6,'r')

	lines1 = infile1.readlines()	
	infile1.close()
	lines2 = infile2.readlines()	
	infile2.close()
	lines3 = infile3.readlines()	
	infile3.close()
	lines4 = infile4.readlines()	
	infile4.close()
	lines5 = infile5.readlines()	
	infile5.close()
	lines6 = infile6.readlines()	
	infile6.close()

	data1  = []
	data2  = []
	data3  = []
	data4  = []
	data5  = []
	data6  = []

	for line in lines1:
		if line.startswith('///') or len(line) < 4 :
			continue
		data1.append(list(map(float,line.split())))
	for line in lines2:
		if line.startswith('///') or len(line) < 4 :
			continue
		data2.append(list(map(float,line.split())))
	for line in lines3:
		if line.startswith('///') or len(line) < 4 :
			continue
		data3.append(list(map(float,line.split())))
	for line in lines4:
		if line.startswith('///') or len(line) < 4 :
			continue
		data4.append(list(map(float,line.split())))
	for line in lines5:
		if line.startswith('///') or len(line) < 4 :
			continue
		data5.append(list(map(float,line.split())))
	for line in lines6:
		if line.startswith('///') or len(line) < 4 :
			continue
		data6.append(list(map(float,line.split())))

	# use latex for font rendering
	mpl.rcParams['text.usetex'] = True
	mpl.rcParams['legend.edgecolor'] = 'white'  

	fig, ax = plt.subplots()
	# plt.title(r'$\textrm{EPPS16, CT14Lo}$', size='xx-large')
	#ax.set_xlabel(r'$x$', size='xx-large')
	#ax.set_ylabel(r'$R_g * f_g$', size='xx-large')
	#plt.ylim(0,10)
	plt.xlim(-8,8)
	plt.grid(visible=True, which='both', axis='both')
	ax.plot(y_points, np.transpose(data1), 'b-', label=r'$pp$', linewidth=1)
	ax.plot(y_points, np.transpose(data2), 'r-', label=r'$AA$', linewidth=1)
	ax.plot(y_points, np.transpose(data3), 'g-', label=r'$sAA$', linewidth=1)
	ax.plot(y_points, np.transpose(data4), 'b--', label=r'$pp\;MC$', linewidth=1)
	ax.plot(y_points, np.transpose(data5), 'r--', label=r'$AA\;MC$', linewidth=1)
	ax.plot(y_points, np.transpose(data6), 'g--', label=r'$sAA\;MC$', linewidth=1)
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
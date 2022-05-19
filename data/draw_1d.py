import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

y_walls = np.linspace(-7.25, 7.25, 30)

obs_names = ['dETdy','dEdy','dETdeta','dEdeta']

for name in obs_names:
	infile_name1  = name + '_sim_pp.dat'
	infile_name2  = name + '_sim_AA.dat'
	infile_name3  = name + '_sim_sAA.dat'

	pdf_name  = name+'.pdf'

	infile1  = open(infile_name1,'r')
	infile2  = open(infile_name2,'r')
	infile3  = open(infile_name3,'r')

	lines1 = infile1.readlines()	
	infile1.close()
	lines2 = infile2.readlines()	
	infile2.close()
	lines3 = infile3.readlines()	
	infile3.close()

	data1  = []
	data2  = []
	data3  = []

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
	ax.plot(y_walls, np.transpose(data1), 'b-', label=r'$pp$', linewidth=1)
	ax.plot(y_walls, np.transpose(data2), 'r-', label=r'$AA$', linewidth=1)
	ax.plot(y_walls, np.transpose(data3), 'g-', label=r'$sAA$', linewidth=1)
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

y_walls = np.linspace(0, 20, 20)

obs_names = ['dETdb','dEdb']

for name in obs_names:
	infile_name1  = name + '_sim_pp.dat'
	infile_name2  = name + '_sim_AA.dat'
	infile_name3  = name + '_sim_sAA.dat'

	pdf_name  = name+'.pdf'

	infile1  = open(infile_name1,'r')
	infile2  = open(infile_name2,'r')
	infile3  = open(infile_name3,'r')

	lines1 = infile1.readlines()	
	infile1.close()
	lines2 = infile2.readlines()	
	infile2.close()
	lines3 = infile3.readlines()	
	infile3.close()

	data1  = []
	data2  = []
	data3  = []

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

	# use latex for font rendering
	mpl.rcParams['text.usetex'] = True
	mpl.rcParams['legend.edgecolor'] = 'white'  

	fig, ax = plt.subplots()
	# plt.title(r'$\textrm{EPPS16, CT14Lo}$', size='xx-large')
	#ax.set_xlabel(r'$x$', size='xx-large')
	#ax.set_ylabel(r'$R_g * f_g$', size='xx-large')
	#plt.ylim(0,10)
	plt.xlim(0,20)
	plt.grid(visible=True, which='both', axis='both')
	ax.plot(y_walls, np.transpose(data1), 'b-', label=r'$pp$', linewidth=1)
	ax.plot(y_walls, np.transpose(data2), 'r-', label=r'$AA$', linewidth=1)
	ax.plot(y_walls, np.transpose(data3), 'g-', label=r'$sAA$', linewidth=1)
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
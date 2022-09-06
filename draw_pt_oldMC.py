import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

#y_walls = np.linspace(-7.25, 7.25, 30)
y_binsize = 0.217949 - (-0.217949) # Indexes 19,20
pt_points = [ 0.75, 1.12536, 1.68858, 2.53369, 3.80175, 5.70446, 8.55943, 12.8433, 19.2711, 28.9159, 43.3878, 65.1026, 97.6852, 146.575, 219.933, 330.005, 495.166, 742.988, 1114.84, 1672.8]

obs_name = ['dNdpT_y=0']

for name in obs_name:
	infile_name1  = 'dNdpTdy_sim_pA_100k_mb_ED_fixxd.dat'
	infile_name2  = 'dNdpTdy_sim_pA_100k_mb_MC_ED_fixxd.dat'

	pdf_name  = name+'.pdf'

	infile1  = open(infile_name1,'r')
	infile2  = open(infile_name2,'r')

	lines1 = infile1.readlines()	
	infile1.close()
	lines2 = infile2.readlines()	
	infile2.close()

	data1  = []
	data2  = []

	n1 = float(lines1[0].split()[2])
	n2 = float(lines2[0].split()[2])

	for line in lines1:
		if line.startswith('///') or len(line) < 4 :
			continue
		data1.append(float(line.split()[21]))
	for line in lines2:
		if line.startswith('///') or len(line) < 4 :
			continue
		data2.append(float(line.split()[21]))

	data7 = [ (m*n2)/(d*n1) if d>0 else 1 for (m,d) in zip(data2, data1)]

	# use latex for font rendering 
	mpl.rcParams['text.usetex'] = True
	mpl.rcParams['legend.edgecolor'] = 'white'  

	fig, ax = plt.subplots()
	plt.xlim(2.5,2510)
	ax.set_xscale('log')
	plt.grid(visible=True, which='both', axis='both')
	ax.plot(pt_points, np.transpose(data1), 'b-', label=r'$pp$', linewidth=1)
	ax.plot(pt_points, np.transpose(data2), 'b--', label=r'$pp\;MC$', linewidth=1)
	ax.minorticks_on()
	ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
	plt.legend(loc='upper right', ncol=1, fontsize='x-large')
	pp = PdfPages(pdf_name)
	plt.savefig(pp, format='pdf')
	pp.close()

	print('Printed to ' + pdf_name)
	plt.show()

	fig, ax = plt.subplots()
	plt.xlim(2.5,2510)
	ax.set_xscale('log')
	plt.grid(visible=True, which='both', axis='both')
	ax.plot(pt_points, np.transpose(data7), 'b-', label=r'$pp$', linewidth=1)
	ax.minorticks_on()
	ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
	plt.legend(loc='upper right', ncol=1, fontsize='x-large')
	pdf_name = name+'_rel.pdf'
	pp = PdfPages(pdf_name)
	plt.savefig(pp, format='pdf')
	pp.close()

	print('Printed to ' + pdf_name)
	plt.show()

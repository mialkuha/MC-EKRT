import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

#y_walls = np.linspace(-7.25, 7.25, 30)
y_binsize = 0.192756 - (-0.192756) # Indexes 19,20
pt_points = [ 2.72832, 5.3985, 7.59385, 10.682, 15.0259, 21.1363, 29.7316, 41.8222, 58.8296, 82.7532, 116.406, 163.743, 230.33, 323.996, 455.753, 641.089, 901.793, 1268.52, 1784.37, 2510 ]

obs_name = ['dNdpT_y=0']

for name in obs_name:
	infile_name1  = 'sigma1jet_sim_pA_100k_mb_ED.dat'
	infile_name2  = 'sigma1jet_sim_pA_100k_mb_MC_ED.dat'

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

	data3 = [ (m*n2)/(d*n1) if d>0 else 1 for (m,d) in zip(data2, data1)]

	# use latex for font rendering 
	mpl.rcParams['text.usetex'] = True
	mpl.rcParams['legend.edgecolor'] = 'white'  

	fig, ax = plt.subplots()
	plt.xlim(2.5,2510)
	ax.set_xscale('log')
	plt.grid(visible=True, which='both', axis='both')
	ax.plot(pt_points, np.transpose(data1), 'b-', label=r'$No\;cons$', linewidth=1)
	ax.plot(pt_points, np.transpose(data2), 'r-', label=r'$Mom\;cons$', linewidth=1)
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
	ax.plot(pt_points, np.transpose(data3), 'b-', label=r'$pA$', linewidth=1)
	ax.minorticks_on()
	ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
	plt.legend(loc='upper right', ncol=1, fontsize='x-large')
	pdf_name = name+'_rel.pdf'
	pp = PdfPages(pdf_name)
	plt.savefig(pp, format='pdf')
	pp.close()

	print('Printed to ' + pdf_name)
	plt.show()
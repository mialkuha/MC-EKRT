import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

#y_walls = np.linspace(-7.25, 7.25, 30)
y_binsize = 0.192756 - (-0.192756) # Indexes 19,20
pt_points = [ 2.72832, 3.46872, 4.41004, 5.60682, 7.12837, 9.06283, 11.5222, 14.6491, 18.6245, 23.6787, 30.1045, 38.2741, 48.6608, 61.8661, 78.655, 100, 200, 376.436, 708.52, 1333.56 ]
#pt_points = [ 2.72832, 5.3985, 7.59385, 10.682, 15.0259, 21.1363, 29.7316, 41.8222, 58.8296, 82.7532, 116.406, 163.743, 230.33, 323.996, 455.753, 641.089, 901.793, 1268.52, 1784.37, 2510 ]

obs_name = ['dNdpT_y=0_2.5m']

for name in obs_name:
	infile_names = [ 'sigma1jet_sim_pA_2500k_mb_PDF.dat',
			 	     'sigma1jet_sim_pA_2500k_mb_nPDF.dat',
			 	     'sigma1jet_sim_pA_2500k_mb_snPDF.dat',
			 	     'sigma1jet_sim_pA_2500k_mb_PDF_MC.dat',
			 	     'sigma1jet_sim_pA_2500k_mb_nPDF_MC.dat',
			 	     'sigma1jet_sim_pA_2500k_mb_snPDF_MC.dat',
			 	     'sigma1jet_sim_pA_2500k_mb_PDF_MC_ND.dat',
			 	     'sigma1jet_sim_pA_2500k_mb_nPDF_MC_ND.dat',
			 	     'sigma1jet_sim_pA_2500k_mb_snPDF_MC_ND.dat' ]

	pdf_name  = name+'.pdf'

	infiles = [ open(n,'r') for n in infile_names ]

	liness = [ f.readlines() for f in infiles ]

	for f in infiles:
		f.close()

	total_ns = [ float(lines[0].split()[2]) for lines in liness ]

	datas = [[float(line.split()[21]) for line in lines if (not line.startswith('///') and len(line) > 3)] for lines in liness]

	data_rel_0 = [ (m*total_ns[3])/(d*total_ns[0]) if d>0 else 1 for (m,d) in zip(datas[3], datas[0])]
	data_rel_1 = [ (m*total_ns[4])/(d*total_ns[1]) if d>0 else 1 for (m,d) in zip(datas[4], datas[1])]
	data_rel_2 = [ (m*total_ns[5])/(d*total_ns[2]) if d>0 else 1 for (m,d) in zip(datas[5], datas[2])]
	data_rel_3 = [ (m*total_ns[6])/(d*total_ns[0]) if d>0 else 1 for (m,d) in zip(datas[6], datas[0])]
	data_rel_4 = [ (m*total_ns[7])/(d*total_ns[1]) if d>0 else 1 for (m,d) in zip(datas[7], datas[1])]
	data_rel_5 = [ (m*total_ns[8])/(d*total_ns[2]) if d>0 else 1 for (m,d) in zip(datas[8], datas[2])]

	# use latex for font rendering 
	mpl.rcParams['text.usetex'] = True
	mpl.rcParams['legend.edgecolor'] = 'white'  

	fig, ax = plt.subplots()
	plt.xlim(2.5,200)
	plt.ylim(0.000005,20)
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.grid(visible=True, which='both', axis='both')
	ax.plot(pt_points, np.transpose(datas[0]), 'b-', label=r'$PDF$', linewidth=1)
	ax.plot(pt_points, np.transpose(datas[1]), 'r-', label=r'$nPDF$', linewidth=1)
	ax.plot(pt_points, np.transpose(datas[2]), 'g-', label=r'$snPDF$', linewidth=1)
	ax.plot(pt_points, np.transpose(datas[3]), 'b--', label=r'$PDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(datas[4]), 'r--', label=r'$nPDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(datas[5]), 'g-.', label=r'$snPDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(datas[6]), 'b-.', label=r'$PDF\;MC\;nd$', linewidth=1)
	ax.plot(pt_points, np.transpose(datas[7]), 'r-.', label=r'$nPDF\;MC\;nd$', linewidth=1)
	ax.plot(pt_points, np.transpose(datas[8]), 'g-.', label=r'$snPDF\;MC\;nd$', linewidth=1)
	ax.minorticks_on()
	ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
	plt.legend(loc='lower left', ncol=1, fontsize='x-large')
	pp = PdfPages(pdf_name)
	plt.savefig(pp, format='pdf')
	pp.close()

	print('Printed to ' + pdf_name)
	plt.show()

	fig, ax = plt.subplots()
	plt.xlim(2.5,100)
	plt.ylim(0.925,1.05)
	ax.set_xscale('log')
	plt.grid(visible=True, which='both', axis='both')
	ax.plot(pt_points, np.transpose(data_rel_0), 'b-', label=r'$PDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(data_rel_1), 'r-', label=r'$nPDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(data_rel_2), 'g-', label=r'$snPDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(data_rel_3), 'b--', label=r'$PDF\;MC\;nd$', linewidth=1)
	ax.plot(pt_points, np.transpose(data_rel_4), 'r--', label=r'$nPDF\;MC\;nd$', linewidth=1)
	ax.plot(pt_points, np.transpose(data_rel_5), 'g--', label=r'$snPDF\;MC\;nd$', linewidth=1)
	ax.minorticks_on()
	ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
	plt.legend(loc='upper right', ncol=1, fontsize='x-large')
	pdf_name = name+'_rel.pdf'
	pp = PdfPages(pdf_name)
	plt.savefig(pp, format='pdf')
	pp.close()

	print('Printed to ' + pdf_name)
	plt.show()

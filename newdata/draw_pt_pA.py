import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

#y_walls = np.linspace(-7.25, 7.25, 30)
y_binsize = 0.192756 - (-0.192756) # Indexes 19,20
pt_points = [ 2.72832, 3.46872, 4.41004, 5.60682, 7.12837, 9.06283, 11.5222, 14.6491, 18.6245, 23.6787, 30.1045, 38.2741, 48.6608, 61.8661, 78.655, 100, 200, 376.436, 708.52, 1333.56 ]
#pt_points = [ 2.72832, 5.3985, 7.59385, 10.682, 15.0259, 21.1363, 29.7316, 41.8222, 58.8296, 82.7532, 116.406, 163.743, 230.33, 323.996, 455.753, 641.089, 901.793, 1268.52, 1784.37, 2510 ]

obs_name = ['dNdpT_2.5m']

for name in obs_name:
	infile_name1  = 'sigma1jet_sim_pA_2500k_mb_PDF.dat'
	infile_name2  = 'sigma1jet_sim_pA_2500k_mb_nPDF.dat'
	infile_name3  = 'sigma1jet_sim_pA_2500k_mb_snPDF.dat'
	infile_name4  = 'sigma1jet_sim_pA_2500k_mb_PDF_MC.dat'
	infile_name5  = 'sigma1jet_sim_pA_2500k_mb_nPDF_MC.dat'
	infile_name6  = 'sigma1jet_sim_pA_2500k_mb_snPDF_MC.dat'
	infile_name7  = 'sigma1jet_sim_pA_2500k_mb_PDF_MC_ND.dat'
	infile_name8  = 'sigma1jet_sim_pA_2500k_mb_nPDF_MC_ND.dat'
	infile_name9  = 'sigma1jet_sim_pA_2500k_mb_snPDF_MC_ND.dat'

	pdf_name  = name+'.pdf'

	infile1  = open(infile_name1,'r')
	infile2  = open(infile_name2,'r')
	infile3  = open(infile_name3,'r')
	infile4  = open(infile_name4,'r')
	infile5  = open(infile_name5,'r')
	infile6  = open(infile_name6,'r')
	infile7  = open(infile_name7,'r')
	infile8  = open(infile_name8,'r')
	infile9  = open(infile_name9,'r')

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
	lines7 = infile7.readlines()	
	infile7.close()
	lines8 = infile8.readlines()	
	infile8.close()
	lines9 = infile9.readlines()	
	infile9.close()

	data1  = []
	data2  = []
	data3  = []
	data4  = []
	data5  = []
	data6  = []
	data7  = []
	data8  = []
	data9  = []

	n1 = float(lines1[0].split()[2])
	n2 = float(lines2[0].split()[2])
	n3 = float(lines3[0].split()[2])
	n4 = float(lines4[0].split()[2])
	n5 = float(lines5[0].split()[2])
	n6 = float(lines6[0].split()[2])
	n7 = float(lines7[0].split()[2])
	n8 = float(lines8[0].split()[2])
	n9 = float(lines9[0].split()[2])

	for line in lines1:
		if line.startswith('///') or len(line) < 4 :
			continue
		data1.append(sum(map(float,line.split()[2:])))
	for line in lines2:
		if line.startswith('///') or len(line) < 4 :
			continue
		data2.append(sum(map(float,line.split()[2:])))
	for line in lines3:
		if line.startswith('///') or len(line) < 4 :
			continue
		data3.append(sum(map(float,line.split()[2:])))
	for line in lines4:
		if line.startswith('///') or len(line) < 4 :
			continue
		data4.append(sum(map(float,line.split()[2:])))
	for line in lines5:
		if line.startswith('///') or len(line) < 4 :
			continue
		data5.append(sum(map(float,line.split()[2:])))
	for line in lines6:
		if line.startswith('///') or len(line) < 4 :
			continue
		data6.append(sum(map(float,line.split()[2:])))
	for line in lines7:
		if line.startswith('///') or len(line) < 4 :
			continue
		data7.append(sum(map(float,line.split()[2:])))
	for line in lines8:
		if line.startswith('///') or len(line) < 4 :
			continue
		data8.append(sum(map(float,line.split()[2:])))
	for line in lines9:
		if line.startswith('///') or len(line) < 4 :
			continue
		data9.append(sum(map(float,line.split()[2:])))

	data10 = [ (m*n4)/(d*n1) if d>0 else 1 for (m,d) in zip(data4, data1)]
	data11 = [ (m*n5)/(d*n2) if d>0 else 1 for (m,d) in zip(data5, data2)]
	data12 = [ (m*n6)/(d*n3) if d>0 else 1 for (m,d) in zip(data6, data3)]
	data13 = [ (m*n7)/(d*n1) if d>0 else 1 for (m,d) in zip(data7, data1)]
	data14 = [ (m*n8)/(d*n2) if d>0 else 1 for (m,d) in zip(data8, data2)]
	data15 = [ (m*n9)/(d*n3) if d>0 else 1 for (m,d) in zip(data9, data3)]

	# use latex for font rendering 
	mpl.rcParams['text.usetex'] = True
	mpl.rcParams['legend.edgecolor'] = 'white'  

	fig, ax = plt.subplots()
	plt.xlim(2.5,200)
	plt.ylim(0.000005,20)
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.grid(visible=True, which='both', axis='both')
	ax.plot(pt_points, np.transpose(data1), 'b-', label=r'$PDF$', linewidth=1)
	ax.plot(pt_points, np.transpose(data2), 'r-', label=r'$nPDF$', linewidth=1)
	ax.plot(pt_points, np.transpose(data3), 'g-', label=r'$snPDF$', linewidth=1)
	ax.plot(pt_points, np.transpose(data4), 'b--', label=r'$PDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(data5), 'r--', label=r'$nPDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(data6), 'g-.', label=r'$snPDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(data7), 'b-.', label=r'$PDF\;MC\;nd$', linewidth=1)
	ax.plot(pt_points, np.transpose(data8), 'r-.', label=r'$nPDF\;MC\;nd$', linewidth=1)
	ax.plot(pt_points, np.transpose(data9), 'g-.', label=r'$snPDF\;MC\;nd$', linewidth=1)
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
	ax.plot(pt_points, np.transpose(data10), 'b-', label=r'$PDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(data11), 'r-', label=r'$nPDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(data12), 'g-', label=r'$snPDF\;MC$', linewidth=1)
	ax.plot(pt_points, np.transpose(data13), 'b--', label=r'$PDF\;MC\;nd$', linewidth=1)
	ax.plot(pt_points, np.transpose(data14), 'r--', label=r'$nPDF\;MC\;nd$', linewidth=1)
	ax.plot(pt_points, np.transpose(data15), 'g--', label=r'$snPDF\;MC\;nd$', linewidth=1)
	ax.minorticks_on()
	ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
	plt.legend(loc='upper right', ncol=1, fontsize='x-large')
	pdf_name = name+'_rel.pdf'
	pp = PdfPages(pdf_name)
	plt.savefig(pp, format='pdf')
	pp.close()

	print('Printed to ' + pdf_name)
	plt.show()

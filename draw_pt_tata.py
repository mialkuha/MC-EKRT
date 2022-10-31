import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import csv
from matplotlib.backends.backend_pdf import PdfPages


def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, marker='+')

    # now determine nice limits by hand:
    nof_bins = 20

    xstart = 0
    xend = 0.062
    xbinwidth = (xend-xstart) / (nof_bins - 1)
    xbins = np.arange(xstart, xend, xbinwidth)
    
    ystart = 1
    yend = 20
    ybinwidth = (yend-ystart) / (nof_bins - 1)
    ybins = np.arange(ystart, yend, ybinwidth)

    ax_histx.hist(x, bins=xbins)
    ax_histy.hist(y, bins=ybins, orientation='horizontal')

names = [ 'jets',
    	  'jets_MC',
    	  'jets_SAT',
    	  'jets_SAT_MC']

for name in names:
	infile_name = name + '.csv'

	infile = open(infile_name,'r')
	reader = csv.DictReader(infile)
		
	p_tensors = [(float(jet['E']),float(jet['p_x']),float(jet['p_y']),float(jet['p_z'])) for jet in reader]
	#xs = [float(x) for (x,y,z) in xsys]
	#ys = [float(y) for (x,y,z) in xsys]
	#zs = [float(z) for (x,y,z) in xsys]
	#print(xs)
	#print(ys)

	for (e,x,y,z) in p_tensors:
		m = e*e-(x*x+y*y+z*z)
		if (m >= 0.000001):
			print(m)

	##fig = plt.figure()
	##ax = fig.add_subplot(projection='3d')
	##plt.title('N={}'.format(len(xs)), size='xx-large')
	####ax.set_xlabel(r'$x$', size='xx-large')
	####ax.set_ylabel(r'$R_g * f_g$', size='xx-large')
	##plt.ylim(1,20)
	##plt.xlim(0,0.052)
	###plt.grid(visible=True, which='both', axis='both')
	##ax.scatter(xs, ys, zs, marker='+')
	##ax.minorticks_on()
	##plt.show()

#	# Create a Figure, which doesn't have to be square.
#	fig = plt.figure(constrained_layout=True)
#	# Create the main axes, leaving 25% of the figure space at the top and on the
#	# right to position marginals.
#	ax = fig.add_gridspec(top=0.75, right=0.75).subplots()
#	plt.title('N={}'.format(len(xs)), size='xx-large', loc='right')
#	plt.ylim(1,20)
#	plt.xlim(0,0.062)
#	# The main axes' aspect can be fixed.
#	#ax.set(aspect=1)
#	# Create marginal axes, which have 25% of the size of the main axes.  Note that
#	# the inset axes are positioned *outside* (on the right and the top) of the
#	# main axes, by specifying axes coordinates greater than 1.  Axes coordinates
#	# less than 0 would likewise specify positions on the left and the bottom of
#	# the main axes.
#	ax_histx = ax.inset_axes([0, 1.05, 1, 0.25], sharex=ax)
#	ax_histy = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)
#	# Draw the scatter plot and marginals.
#	scatter_hist(xs, ys, ax, ax_histx, ax_histy)
#
#	### Start with a square Figure.
#	##fig = plt.figure(figsize=(6, 6))
#	### Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
#	### the size of the marginal axes and the main axes in both directions.
#	### Also adjust the subplot parameters for a square plot.
#	##gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
#	##                      left=0.1, right=0.9, bottom=0.1, top=0.9,
#	##                      wspace=0.05, hspace=0.05)
#	###plt.title('N={}'.format(len(xs)), size='xx-large', loc='right')
#	### Create the Axes.
#	##ax = fig.add_subplot(gs[1, 0])
#	###ax.ylim(1,20)
#	###ax.xlim(0,0.052)
#	##ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
#	##ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
#	### Draw the scatter plot and marginals.
#	##scatter_hist(xs, ys, ax, ax_histx, ax_histy)
#
#	pp = PdfPages(pdf_name)
#	plt.savefig(pp, format='pdf')
#	pp.close()
#
#	print('Printed to ' + pdf_name)
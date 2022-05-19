import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import axes3d


y_walls = np.linspace(-7.25, 7.25, 30)
kt_max = 11

#infile_name  = 'sigma1jet_analytical_AA.dat'
#infile_name  = 'sigma1jet_analytical_pp.dat'
infile_name1  = 'sigma1jet_sim_AA.dat'
#infile_name  = 'sigma1jet_sim_pp.dat'
infile_name2  = 'sigma1jet_sim_sAA.dat'
#infile_name  = 'sigmadijet_analytical_AA.dat'
#infile_name  = 'sigmadijet_analytical_pp.dat'
#infile_name  = 'sigmadijet_sim_AA.dat'
#infile_name  = 'sigmadijet_sim_pp.dat'
#infile_name  = 'sigmadijet_sim_sAA.dat'

pdf_name  = infile_name1[:-4]+'_vs_'+infile_name2[14:-3]+'pdf'

infile1  = open(infile_name1,'r')
infile2  = open(infile_name2,'r')

lines1 = infile1.readlines()	
infile1.close()
lines2 = infile2.readlines()	
infile2.close()

data1  = []
data2  = []

for line in lines1:
	if line.startswith('///') or len(line) < 4 :
		continue
	data1.append(tuple(map(float,line.split())))

for line in lines2:
	if line.startswith('///') or len(line) < 4 :
		continue
	data2.append(tuple(map(float,line.split())))

kts = []

for tupl in data1:
	kts.append(tupl[0])
	kts.append(tupl[1])

kts = list(dict.fromkeys(kts))

kts = [kt for kt in kts if kt <= kt_max]

kt_dim = len(kts)
y_dim = len(y_walls)

sigmas1  = []
sigmas2  = []

for tupl in data1:
	sigmas1.append(list(tupl[2:]))
for tupl in data2:
	sigmas2.append(list(tupl[2:]))

sigmas = np.array(sigmas1[:kt_dim]) - np.array(sigmas2[:kt_dim])

# use latex for font rendering
mpl.rcParams['text.usetex'] = True


mpl.rcParams['legend.edgecolor'] = 'white'  

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Grab some test data.

#X, Y, Z = axes3d.get_test_data(0.05)
#X, Y, Z = kts[:-1], y_walls, np.array(sigmas)
dummy = []
for i in range(y_dim):
	dummy.append(kts)
X = np.transpose(dummy)
#print(np.shape(X))
#print(X)
dummy = []
for i in range(kt_dim):
	dummy.append(y_walls)
Y = np.array(dummy)
Z = sigmas
#print(np.shape(Y))
#print(np.shape(Z))
#print(Y)
#print(Z)

# Plot a basic wireframe.
#ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# plt.title(r'$\textrm{EPPS16, CT14Lo}$', size='xx-large')
#ax.set_xlabel(r'$k_T$', size='xx-large')
#ax.set_ylabel(r'$y$', size='xx-large')
plt.ylim(-8,8)
plt.xlim(2.5,11)
ax.set_zlim(-1,4)
#plt.xscale('log')
#plt.grid(visible=True, which='both', axis='both')
#print(np.shape(kts[:-1]))
#print(np.shape(y_walls))
#print(np.shape(np.array(sigmas)))
surf = ax.plot_surface(X, Y, Z, cmap=mpl.cm.coolwarm, linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
#ax.minorticks_on()
#ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)

pp = PdfPages(pdf_name)
plt.savefig(pp, format='pdf')
pp.close()

print('Printed to ' + pdf_name)
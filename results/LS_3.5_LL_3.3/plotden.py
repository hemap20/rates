#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import os
from matplotlib import rcParams 
import matplotlib as mpl 

#rc('text', fontsize=32, usetex=True)
#rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#rc('text', fontsize=28)
#rc('font', weight='bold')



fz1 = 30 
fz2 = 32 
fz3 = 30 
colors = ["purple", "deeppink", "peru", "teal" , "lawngreen", "springgreen", "tomato", "yellow" , "orange"]

#rcParams.update({'figure.autolayout': True}) 
#mpl.rc('axes',linewidth=2) 

plt.figure(figsize=(10, 8))

##################################################

f2 = open("./density_5.txt", 'r')
lines = f2.readlines()
f2.close()

x = []
density = []

for line in lines:
	p=line.split()
	x.append(float(p[0]))
	density.append(float(p[2]))

xv = np.array(x)
densityv = np.array(density)



plt.plot(xv,densityv,'peru',lw=2.0, label=r'$\sigma_{LL} = 3.2 \AA, \epsilon_{LL}= 3.4 kJ/mol,\sigma_{LS} = 3.3 \AA, \epsilon_{LS}= 4$ kJ/mol')



plt.xlabel(r'Distance from bottom ($\AA$)')
plt.ylabel('Density (mol/cm3)')

plt.xlim(0, 70)
plt.ylim(0,0.1)
leg = plt.legend(frameon=False) 
plt.savefig('density_run_5.png', format='png',dpi=300, bbox_inches='tight')
plt.show()
##################################################

# f2 = open("./llepsilon_2.6_ls_3.5/density.txt", 'r')
# lines = f2.readlines()
# f2.close()

# x = []
# density = []

# for line in lines:
# 	p=line.split()
# 	x.append(float(p[0]))
# 	density.append(float(p[2]))

# xv = np.array(x)
# densityv = np.array(density)



# plt.plot(xv,densityv,'deeppink',lw=2.0, label=r'$\sigma_{LL} = 3.2 \AA, \epsilon_{LL}= 2.6 kJ/mol,\sigma_{LS} = 3.3 \AA, \epsilon_{LS}= 3.5$ kJ/mol')


# #################################################


# f2 = open("./llepsilon_2.6_ls_3/density.txt", 'r')
# lines = f2.readlines()
# f2.close()

# x = []
# density = []

# for line in lines:
# 	p=line.split()
# 	x.append(float(p[0]))
# 	density.append(float(p[2]))

# xv = np.array(x)
# densityv = np.array(density)



# plt.plot(xv,densityv,'yellow',lw=2.0, label=r'$\sigma_{LL} = 3.2 \AA, \epsilon_{LL}= 2.6 kJ/mol,\sigma_{LS} = 3.3 \AA, \epsilon_{LS}= 3$ kJ/mol')


# #################################################


# f2 = open("./llepsilon_2.6_ls_2.5/density.txt", 'r')
# lines = f2.readlines()
# f2.close()

# x = []
# density = []

# for line in lines:
# 	p=line.split()
# 	x.append(float(p[0]))
# 	density.append(float(p[2]))

# xv = np.array(x)
# densityv = np.array(density)



# plt.plot(xv,densityv,'teal',lw=2.0, label=r'$\sigma_{LL} = 3.2 \AA, \epsilon_{LL}= 2.6 kJ/mol,\sigma_{LS} = 3.3 \AA, \epsilon_{LS}= 2.5$ kJ/mol')


#################################################





#os.system('display freeEnergy.png')

# (projected) band structure plotting code for VASP 6.3.1:
# --------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

#TODO: total per ion
#TODO: fix dk in plots
#TODO: make plotting, labeling, titling, savefig, ... automatic! CHOICE: ions, orbitals, ...
#TODO: clean code
#TODO: pDOS (DOSCAR from Band calc NEDOS)
#TODO: ISPIN, spinors, ...
#TODO: make input params (N_b, ...) being read from inputfiles
#TODO: params consistent w/ hopping code

#################################################################
E_F = 8.9499947474 # From LWL!

N_bands = 120
N_kseg = 50
K_points = ['$\Gamma$', 'X', 'M', '$\Gamma$']

project = True
N_orbs = 16 # 1, 4, 9, 16, ...

ions = [3, 4]
orbs = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2y2"]

#################################################################

WB = open('input/OUTCAR_PRIM', "r")
WB = WB.read()
WB = np.array([x for x in WB.split(" ") if x != '' and x != '\n'])
k_start = np.where(WB == 'occupation')[0] + 2

E = np.zeros((len(k_start), N_bands))
sum_dk = 0
dk = [0]

for i in range(1,len(k_start)):
    dk += [np.sqrt((float(WB[k_start[i]-9])-float(WB[k_start[i-1]-9]))**2+(float(WB[k_start[i]-8])-float(WB[k_start[i-1]-8]))**2+(float(WB[k_start[i]-7])-float(WB[k_start[i-1]-7]))**2)]
    sum_dk += dk[i]
    e = np.zeros(N_bands)
    for j in range(N_bands):
        e[j] = WB[k_start[i]+3*j]
    E[i] = e
    
dk = (1/sum_dk)*np.array(dk)
E = np.transpose(E)

kk = np.zeros(len(k_start))
for i in range(1,len(k_start)):
    if i == 0:
        kk[0] = 0
    else:
        kk[i] = kk[i-1] + dk[i]
    

K_values = [kk[0]]
for K in range(1, len(K_points)):
    K_values += [kk[K*N_kseg-1]]

E -= E_F

E = E[:,1:] #TODO WHY?
kk = kk[1:] #TODO WHY?

PC = open('input/PROCAR_PRIM', "r")
PC = PC.read()
PC = np.array([x for x in PC.split(" ") if x != '' and x != '\n'])

N_k = int(PC[5])
N_b = int(PC[9])
N_i = int(PC[13])

if project: 
    P = PC[:]

    Data = np.zeros((N_b, N_k, N_i, 3))

    for k in range(N_k):
        for b in range(N_b):
            for i in range(N_i):
                if N_orbs >= 1: Data[b, k, i, 0] = float(P[32 + N_orbs + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k])
                for o in [0, 1, 2]:
                    if N_orbs >= 4: Data[b, k, i, 1] += float(P[33 + N_orbs + o + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k])
                for o in [2]:
                    if N_orbs >= 9: Data[b, k, i, 2] += float(P[36 + N_orbs + o + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k])


# Plotting lm-decomposed band structures:
# ---------------------------------------

fig, axs = plt.subplots()
norm = plt.Normalize(0, 1)
for i in range(len(E)):
    seg = np.array([[[kk[j], E[i][j]], [kk[j+1], E[i][j+1]]] for j in range(len(kk)-1)])

    if project:  p = Data[i, :, 3, 2] + Data[i, :, 4, 2]
    
    lc = LineCollection(seg, cmap='BuPu',norm=norm)
    if project:  lc.set_array(p)
    lc.set_linewidth(1.5)
    line = axs.add_collection(lc)

if project: fig.colorbar(line, ax=axs)

plt.ylabel(r'$E-E_{F}$ [eV]', fontsize=15)
plt.xlabel('K-points', fontsize=15)
plt.ylim(-2,3)
plt.axhline(0, linestyle=(0, (5,5)), linewidth=1, color='g', alpha=1)
plt.xticks(K_values, K_points, fontsize=15)
plt.yticks(fontsize=15)
plt.grid()
# plt.title('Cu: d$_{z^2}$', fontsize=15) #CHANGE THIS
plt.savefig("figures/La3Ni2O7_PRIM_PBE_bands_Ni_d_z2.png")

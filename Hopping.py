import numpy as np
import spglib

np.set_printoptions(threshold=np.inf)

#######################################################
path_to_POSCAR = 'input/POSCAR_PRIM'
path_to_KPOINTS = 'input/KPOINTS_PRIM'
path_to_OUTCAR = 'input/OUTCAR_LWL_PRIM'
path_to_WANPROJ = 'input/WANPROJ_PRIM'

N_b = 120 #TODO
N_w = 2 #TODO

nn = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]]) # Translation R along POSCAR-defined a,b,c!

#TODO: clean all the bs code below... make some plotting display of hopping matrices in nice gui?

#######################################################

def read_poscar(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    scaling_factor = float(lines[1].strip())

    lattice = []
    for i in range(2, 5):
        lattice.append([float(x) for x in lines[i].strip().split()])
    lattice = np.array(lattice) * scaling_factor

    element_counts = [int(x) for x in lines[6].strip().split()]

    numbers = []
    atom_type = 1
    for count in element_counts:
        numbers.extend([atom_type] * count)
        atom_type += 1

    fractional_positions = []
    for line in lines[8:]:
        fractional_positions.append([float(x) for x in line.strip().split()])

    return lattice, fractional_positions, numbers

def read_kpoints(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    mesh = [int(x) for x in lines[3].strip().split()]
    
    return mesh

cell = read_poscar(path_to_POSCAR)
print("Space group:", spglib.get_spacegroup(cell, symprec=1e-5))

kmesh = read_kpoints('input/KPOINTS_PRIM')

N_k = kmesh[0]*kmesh[2]*kmesh[2]
print("Number of reducible kpoints:", N_k)

mapping, grid = spglib.get_ir_reciprocal_mesh(kmesh, cell, is_shift=[0, 0, 0])

equiv = np.zeros((N_k, 2, 3))
for i, (ir_gp_id, gp) in enumerate(zip(mapping, grid)):
    # print("%3d ->%3d %s -> %s" % (i, ir_gp_id, gp.astype(float) / kmesh, grid[ir_gp_id].astype(float) / kmesh))
    equiv[i, 0] = gp.astype(float) / kmesh
    equiv[i, 1] =  grid[ir_gp_id].astype(float) / kmesh

N_r = len(np.unique(mapping))
print("Number of irreducible kpoints: %d" % N_r)

WB = open(path_to_OUTCAR, "r")
WB = WB.read()
WB = np.array([x for x in WB.split(" ") if x != '' and x != '\n'])
k_start = np.where(WB == 'occupation')[0]+2

ind = np.where(WB == 'Fermi')[0][-1]
E_F = float(WB[ind+2])

k_points= np.zeros((N_r, 3))
energies = np.zeros((N_r, N_b))
for i in range(N_r):
    k = [float(WB[k_start[i]-9]), float(WB[k_start[i]-8]), float(WB[k_start[i]-7][:-1])]
    k_points[i] = k
    for j in range(N_b):
        Band = j
        energies[i, j] = WB[k_start[i]+3*j]

PROJ = open(path_to_WANPROJ, "r")
PROJ = PROJ.read()
PROJ = np.array([x for x in PROJ.split(" ") if x != '' and x != '\n'])[15+N_k*4:]

PROJ2 = open(path_to_WANPROJ, "r")
PROJ2 = PROJ2.read()
PROJ2 = np.array([x for x in PROJ2.split(" ") if x != '' and x != '\n'])[15:]

KPOINTS = np.zeros((N_k, 3))
for i in range(N_k):
    KPOINTS[i] = np.array([round(float(PROJ2[1+4*i]), 4), round(float(PROJ2[2+4*i]), 4), round(float(PROJ2[3+4*i]), 4)])

KP = np.zeros((N_k, 3))
for i in range(N_k):
    KP[i] = np.array([float(PROJ2[1+4*i]), float(PROJ2[2+4*i]), float(PROJ2[3+4*i])])
    o = np.where(abs(KPOINTS[i]) > 0.5)[0]
    KP[i][o] = -1*np.sign(KPOINTS[i][o])*(1.0-abs(KPOINTS[i][o]))
    KP[i] = np.round(KP[i], 4)

KPOINTS2 = np.zeros((N_k, 3))
for k in range(len(equiv)):
    for i in range(3):
        equiv[k,0,i] = round(equiv[k,0,i], 4)
        equiv[k,1,i] = round(equiv[k,1,i], 4)
        KPOINTS2[k,i] = round(KP[k,i],4)


ENERGIES = np.zeros((N_k, N_b))
for k in range(len(KPOINTS)):

    iii = np.where(np.all(equiv[:,0]==KPOINTS2[k],axis=1))[0]
    ii = np.where(np.all(k_points==equiv[iii,1],axis=1))[0]

    ENERGIES[k] = energies[ii[0]]
    

T = np.zeros((N_k, N_b, N_w),dtype=complex)
for K in range(len(KPOINTS)):
    for b in range(N_b):
        for i in range(N_w):
            Re, Im = float(PROJ[7+4*N_w*b+4*i+(4*N_w*N_b+5)*K]), float(PROJ[8+4*N_w*b+4*i+(4*N_w*N_b+5)*K])
            T[K, b, i] = complex(Re,Im)

t = np.zeros((len(nn), N_w, N_w))
for R in range(len(nn)):    
    for i in range(N_w):
        for j in range(N_w):
            for n in range(N_b):
                for k in range(len(KP)):
                    t[R, i, j] -= ((np.conj(T[k, n, i])*(ENERGIES[k, n]-E_F)*T[k, n, j])*np.exp(1j*2*np.pi*np.dot(KP[k], nn[R]))).real
                
    print('Hopping matrix to', nn[R],':\n \n', 't_ij=', (1/N_k)*t[R,0], '\n')
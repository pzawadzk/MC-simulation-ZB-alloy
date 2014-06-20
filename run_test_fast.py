#!/usr/bin/env python
from random import choice
from zb import czts_kesterite
from acc_mc_disorder import run_mc
import numpy as np
import os
from pylada.crystal import supercell

dir0 = 'test'

Niter_melt = 1000

######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########
# Hamiltonian
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########

dh = np.array([-0.32752784, -0.29256123, -0.59147627, -0.87821934, -1.08734453,
       -0.26111739, -0.55886433, -0.78891833, -0.90441077, -0.46754391,
       -0.66754745, -0.67916151, -0.48116762, -0.40339855, -0.23411278])
#acccure czts phase transtion
dh = np.array([-0.31233548, -0.26810176, -0.579079  , -0.8969619 , -1.086384  ,
       -0.22762914, -0.54743073, -0.77016516, -0.89829263, -0.45909058,
       -0.6633367 , -0.66669739, -0.5115215 , -0.41838264, -0.20690877])

######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########
#CELL
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########
#[60.000000000000007, 60.000000000000007, 60.000000000000007], [10.606601717798213, 10.606601717798213, 10.606601717798213]]
str = czts_kesterite()

cell = str.cell
det = np.linalg.det(cell)
Ts = [-200.]*24
N = 2
scell = np.array([[2*N, 0, 0], [0, 2*N, 0], [0, 0, N]])
#scell = np.array([[14, 0, 0], [0, 14, 0], [0, 0, 7]])
cell = np.dot(cell, scell)
print cell

structure  = supercell(str, cell)

types = ['Cu', 'Zn', 'Sn']
sites = {}


for type in types:
  sites[type] = []

det_tot = 0
for id, atom in enumerate(structure):
  if atom.type in types:
    sites[atom.type].append(id)
    det_tot += 1


sites0 = [ sites['Cu'], sites['Zn'], sites['Sn'] ] 
N0 = [ len(sites[type]) for type in types]

print 'DET', det_tot
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########
#for N_Sn in range(det+1): 
#N_Sn takes valus from 0-10
lattices = []
if not os.path.exists(dir0):
  os.makedirs(dir0)

os.chdir(dir0)
Tm = 1500 
T =  -200
seed = 0
if 1:

    N_Cu = det_tot//2
    N_Sn = det_tot//4
    N_Zn = det_tot - N_Sn - N_Cu
    sites = sites0[:]

    N = [N_Cu, N_Zn, N_Sn]
    dN = [N[i]-N0[i] for i in range(3)]

    excess = []
    for ii in range(3):
      if dN[ii] <0:
        for ni in range(abs(dN[ii])):
          id = choice(sites[ii])
          id_index = sites[ii].index(id)
          sites[ii].pop(id_index)
          excess.append(id)

    for ii in range(3):
      if dN[ii] >0:
        for ni in range(dN[ii]):
          id = choice(excess)
          id_index = excess.index(id)
          excess.pop(id_index)
          sites[ii].append(id)

    for ii in range(3):
      for ni in sites[ii]:
        structure[ni].type = types[ii]

    path = 'cu_sn_%d_%d_%d_%d_T_%.1f_seed_%d_mc_4_fit_large_cell'%(N_Cu, N_Sn, N_Zn, det_tot, T, seed)
    if not os.path.exists(path):
      os.mkdir(path)
    os.chdir(path)
    niter = int((Tm-T)/100)+1
    restart = False
    T = np.linspace(Tm, T, niter) +273.
    T = None #np.linspace(Tm, T, niter) +273.
    print T
    latt = run_mc(structure, dh, types, T, Niter_melt, tag='', restart = restart, seed = seed)
    lattices.append(latt)

    os.chdir('..')
os.chdir('..')

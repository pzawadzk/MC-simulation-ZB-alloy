#!/usr/bin/env python
import time
from sets import Set
from pylada.crystal import write, read, supercell, neighbors
import numpy as np
from random import choice
import random

def get_types_in_1st_cn(neight_list, all_atoms, types_index, all_atoms_tmp):
  """ Calculates composition of the coordination shell. """
  atoms = [0, 0, 0]
  for i in neight_list:
    if i in all_atoms_tmp.keys(): 
      symb = all_atoms_tmp[i]
    else:
      symb = all_atoms[i]
    atoms[types_index[symb]] += 1
  return atoms

def calc_motif(neight_list, all_atoms, types_index, all_atoms_tmp={}):
  """ Calculates type of motif assisted with a particular composition of the coordination shell. """
  matrix  = [ [0,1,2,3,4], [5,6,7,8,-1], [9,10,11,-1,-1], [12,13,-1,-1,-1], [14,-1,-1,-1,-1],]
  atoms = get_types_in_1st_cn(neight_list, all_atoms, types_index, all_atoms_tmp)
  el = matrix[atoms[0]][atoms[1]]
  return el 

def check_composition(dn_tot):
  """ Calculated the total composition from total number of motifs. """
  p = np.zeros([5,5])
  p[0]= dn_tot[0:5]
  p[1,:4]= dn_tot[5:9]
  p[2,:3]= dn_tot[9:12]
  p[3,:2]= dn_tot[12:14]
  p[4,:1]= dn_tot[14]

  N_A, N_B, N_C = 0, 0, 0
  for i in range(5):
    for j in range(0,5-i):
      N_A += i*p[i,j]/4.
      N_B += j*p[i,j]/4.
      N_C += (4-i-j)*p[i,j]/4.
  return N_A, N_B, N_C

def run_mc(structure, dh, types, T = None, Niter = 1000, tag='', restart= False, seed = 0, output_freq=100, output_freq_dat=100):
  """ Runs Monte Carlo Simulations.

      Args:
        structure (pyjama structure object): structure  
        dh (list): Model Hamiltonian
        T0 (np.array): Temperature in Kelvin
        Niter (int): Number of iterations per T step
        tag (str): Name of the output file .dat and .start
        restart (bool): Restart from existing .vasp and .dat file
        seed (int): Random seed
        output_freq (int): Frequency with which structures are written to .vasp
        output_freq_dat (int): Frequency with which data is written to .dat

  """
  Niter = int(Niter)
  random.seed(seed)
  np.random.seed(seed)
  if T == None:
    T = np.array([2000])
  beta = 1./(8.617e-5 * T)
  Nstep = len(T)

  tag +=  'N_%s_%s_%s_'%(types[0], types[1], types[2])
  types_index = {types[0]:0, types[1]:1, types[2]:2}

  file_dat = 'mc_%s.dat'%tag
  if restart:
    try:
      fdat = open(file_dat, 'r')
      mc_f = np.loadtxt(fdat)
      last_iter = int(mc_f[-1][0])
      fdat.close()
      structure = read.poscar('poscar_%010d.vasp'%last_iter)
      last_iter += 1
    except:
      print 'ERROR restarting from file:', file_dat
  else:
    last_iter = 1
    fdat = open(file_dat, 'w')
    fdat.close()

  #list of all sites that can be exchanged in MC
  all_MC_sites = []
  # list of all atoms
  all_atoms = {}
  all_index = []
  for i, atom in enumerate(structure):
    if atom.type in types:
      all_MC_sites.append(i)
    #all_atoms.append(atom.type)
    all_atoms[i] = atom.type
    all_index.append(atom)

  ########################
  start_time = time.time()
  # dictionarry of list of 1st neighbours
  atoms_dict = {}
  S_motifs = {}
  for i, atom in enumerate(structure):
    atoms_dict[i] = []

    neighs = [n for n in neighbors(structure, 4, atom.pos)]
    for n in neighs:
      atoms_dict[i].append(all_index.index(n[0]))

    if atom.type == 'S':
      S_motifs[i] = calc_motif(atoms_dict[i], all_atoms, types_index)
  end_time = time.time()
  print 'Initialization of distance matrix'
  print("Elapsed time was %g seconds" % (end_time - start_time))
  ########################

  # initial number of motifs
  dn_tot = np.zeros(15)
  for index in S_motifs.keys():
    dn_tot[S_motifs[index]] += 1

  param_size = len(dn_tot)
  # initial energy
  dE_tot =  2*np.dot(dn_tot, dh)

  # initial energy
  f_start= open('mc_%s.start'%tag, 'a')
  print >>f_start, 'Start Motives:'
  out_to_file = ''
  for k in range(param_size):
      out_to_file  += '\t%.3f' %dn_tot[k]
  print >>f_start, out_to_file
  print >>f_start, check_composition(dn_tot)
  f_start.close()

  ########################
  start_time = time.time()
  for iter in xrange(last_iter, Niter*Nstep+last_iter):

    id_1 = choice(all_MC_sites)
    id_2 = choice(all_MC_sites)
    while all_atoms[id_1] == all_atoms[id_2]:
      id_1 = choice(all_MC_sites)
      id_2 = choice(all_MC_sites)

    all_atoms_tmp = {}
    all_atoms_tmp[id_1] = all_atoms[id_2]
    all_atoms_tmp[id_2] = all_atoms[id_1]

    S_to_update = Set(atoms_dict[id_1] + atoms_dict[id_2])

    dn = np.zeros(15)
    S_motifs_tmp = {}
    for S_index in S_to_update:
      dn[S_motifs[S_index]] -= 1
      S_motifs_tmp[S_index] = calc_motif(atoms_dict[S_index], all_atoms, types_index, all_atoms_tmp)
      dn[S_motifs_tmp[S_index]] += 1
    dE =  2*np.dot(dn, dh)
    end_time = time.time()

    if dE<0:
      accept = 1
    else:
      A_move = np.exp( - beta[(iter-last_iter)//Niter] * dE)
      if A_move > np.random.random():
        accept = A_move
      else:
        accept = 0
    if accept:
      dE_tot += dE
      dn_tot += dn
      for S_index in S_to_update:
        S_motifs[S_index] = S_motifs_tmp[S_index]
      for id in [id_1, id_2]:
        all_atoms[id] = all_atoms_tmp[id]

    if iter % output_freq_dat == 0:
      fdat = open(file_dat, 'a')
      out_to_file  = '%d'%iter
      for k in range(param_size):
        out_to_file  += '\t%.3f' %dn_tot[k] 
      out_to_file  += '\t%.3f\t%.3f\t%.3f\t%d\t%.3f'% (dE_tot, dE, accept,  iter, T[(iter-last_iter)//Niter])
      print >> fdat, out_to_file
      fdat.close()

    if iter % output_freq == 0:
      for j, atom in all_atoms.items():
      #for j, atom in enumerate(all_atoms):
        structure[j].type = atom
      write.poscar(structure, 'poscar_%010d.vasp'%iter, vasp5=True)

  end_time = time.time()
  print("Write Elapsed time was %g seconds per iter" % ((end_time - start_time)/(Niter*Nstep)))
  # initial energy
  f_start= open('mc_%s.start'%tag, 'a')
  print >>f_start, 'Final Motives:'
  out_to_file = ''
  for k in range(param_size):
      out_to_file  += '\t%.3f' %dn_tot[k]
  print >>f_start, out_to_file
  print >>f_start, check_composition(dn_tot)
  f_start.close()

  ########################


from pylada.vasp import Extract
import pickle
import numpy as np
import sys
import os
from str_util import  get_index
from scipy.optimize import leastsq

def count_configurations_5(structure, atom_1, atoms, n_neigh=4, index_list=None, if2nd=False):
  from pylada.crystal import neighbors
  conf = {}
  S_2nd = {}
  N = len(atoms)
  conf_count = np.zeros([5, 5])
  if index_list == None:
    index_list = range(len(structure))
  for id in index_list:
    atom = structure[id]
    if atom.type == atom_1:
      S_2nd[id] = []
      p = atom.pos
      neighs = [n for n in neighbors(structure, n_neigh, p)]
      count = np.zeros(N)
      for n in neighs:
        neighs_2nd = [n_2nd for n_2nd in neighbors(structure, n_neigh, n[0].pos)]
        for n_2nd in neighs_2nd:
          S_2nd[id].append(get_index(structure, n_2nd[0]))
        for i, atom in enumerate(atoms):
          if n[0].type == atom:
            count[i]  +=1
        if sum(count) == 4:
          break
          #print n.distance, bulk.atoms[n.index].type
      tag = ''
      for i, atom in enumerate(atoms):
        tag += atom +'%d_'%count[i]
      if tag not in conf.keys():
        conf[tag] = 0 
      conf[tag] +=1
      conf_count[int(count[0]),int(count[1])] +=1
  if if2nd:
    conf_count = np.concatenate([conf_count[0], conf_count[1,:-1], conf_count[2,:-2],conf_count[3,:-3],conf_count[4,:-4]])
    return conf, conf_count.tolist(), S_2nd
  else:
    conf_count = np.concatenate([conf_count[0], conf_count[1,:-1], conf_count[2,:-2],conf_count[3,:-3],conf_count[4,:-4]])
    return conf, conf_count.tolist()

def get_Hf(calc):

  dict ={'Cu':[], 'Sn':[], 'Sb':[], 'S':[], 'Zn':[], 'Mg':[], 'P': [], 'Al':[], 'N':[], 'O':[], 'As': [], 'Bi': [], 'Si':[], 'In':[], 'Se':[], 'Ga':[]}

  for id, atom in enumerate(calc.structure):
    if atom.type in dict.keys():
      dict[atom.type].append(id)
    else:
      dict[atom.type].append(id)

  #print dict
  mu = {'Cu': -1.97, 'Sn':  -3.79, 'S':  -4.00 , 'Sb': -4.29, 'Zn':-0.84, 'Mg': -0.99, 'P':-5.64, 'Al': -3.02, 'N':-8.51, 'O':-4.76, 'As': -5.06, 'Bi': -4.39, 'Si':-4.99, 'In': -2.31, 'Se': -3.55, 'Ga': -2.37}

  import numpy as np
  E = float(calc.energy_sigma0)
  N = len(calc.structure)
  Ef = 0
  for key in dict.keys():
    Ef -= mu[key]*len(dict[key])


  Ef += E
  Ef /=N
  return Ef

def extract_motifs(dirs0, center, atoms, calc_mad=False):
  x = []
  E = []
  Formula = []
  for dir0 in dirs0:
    for dir in os.listdir(dir0):
      if os.path.isdir(os.path.join(dir0,dir)):
        #print os.path.join(dir0,dir)
        file_data = os.path.join(dir0, dir, 'data.pckl')
        if 1: #not os.path.isfile(file_data):
            if 'non-magnetic' in os.listdir(os.path.join(dir0, dir)):
                file = os.path.join(dir0, dir, 'non-magnetic', 'OUTCAR')
            elif 'relax' in os.listdir(os.path.join(dir0, dir)):
                file = os.path.join(dir0, dir, 'relax', 'OUTCAR')
            else:
                file = os.path.join(dir0, dir, 'OUTCAR')
            #if os.path.isfile(file):
            try:
                #print file
                calc = Extract(file)
                E0f = get_Hf(calc)
                N = len(calc.structure) 
                if calc_mad:
                  #str_ase = structure_lada_to_ase(calc.structure)
                  #U = get_lattice_energy(str_ase, q_dict)
                  #U = np.sum(U)
                  U = 0
                else:
                  U = 0
                count, conf_count = count_configurations_5(calc.structure, center, atoms = atoms, n_neigh = 8)
                #print count
                conf_count  = np.array(conf_count)
                n = conf_count/np.sum(conf_count)
                f=open(file_data, 'w')
                pickle.dump([E0f, n, U],f)
                f.close()
                E.append(E0f)
                x.append(n)
                #print dir, E0f, U
            #else:
            except:
                #print 'Error', dir
                pass
        else:
            f = open(file_data, 'r')
            [E0f, n, U] = pickle.load(f)
            f.close()

            E.append(E0f)
            x.append(n)
    
  return x, E#, Formula

def fit_motifs(x, E, p0 = None, indexes = None):
#    x = np.array(x)
#    E = np.array(E)
    if indexes == None:
      indexes = range(len(x[0]))

    N = len(indexes)
    if p0 == None: 
      p0 = np.zeros(15) - 0.5

    p = p0[indexes] 
    #print p


    p = leastsq(residuals, p, args=(E, x, p0, indexes ))#, full_output=True)
    for count, index in enumerate(indexes):
      p0[index] = p[0][count]
        
    Efit = func(p0,x, p0, range(15))
    
    #dH= abs(E - Efit)
    #print 'In sample MAE', 1000*np.mean(dH), 'meV'    

    return p0, Efit

def func(p, x, p0, indexes):
  f = 0
  for i in range(15):
    if  i in indexes: 
      f += p[indexes.index(i)] * x[:,i] 
    else:
      f += p0[i] * x[:,i]
  return f

def residuals(p, y, x, p0, indexes): 
  err = y - func(p, x, p0, indexes)
  return err




def print_stable_motifs(dH):
  dh_dict = {}

  p = np.zeros([5,5])
  p[0]= dH[0:5]
  p[1,:4]= dH[5:9]
  p[2,:3]= dH[9:12]
  p[3,:2]= dH[12:14]
  p[4,:1]= dH[14]

  stable = []
  unstable = []
  for i in range(5):
    for j in range(5-i):
      stable.append([i,j])

  for i in range(5):
    for j in range(5-i):
      N_A = i
      N_B = j
      N_C = (4-i-j)
      if N_A>0 and N_A<4 and N_C>0 and N_C<4:
        dh = p[i-1,j] + p[i+1, j] - 2 * p[i,j] 
        #print '2 %d_%d_%d -> %d_%d_%d + %d_%d_%d\t\t%d meV' %(N_A, N_B, N_C, N_A-1, N_B, N_C+1,N_A+1, N_B, N_C-1, 2000*dh) 
      if N_B>0 and N_B<4 and N_C>0 and N_C<4:
        dh = p[i,j-1] + p[i, j+1] - 2 * p[i,j] 
        #print '2 %d_%d_%d -> %d_%d_%d + %d_%d_%d\t\t%d meV' %(N_A, N_B, N_C, N_A, N_B-1, N_C+1,N_A, N_B+1, N_C-1, 2000*dh) 
      if N_A>0 and N_A<4 and N_B>0 and N_B<4:
        dh = p[i+1,j-1] + p[i-1, j+1] - 2 * p[i,j] 
        #print '2 %d_%d_%d -> %d_%d_%d + %d_%d_%d\t\t%d meV' %(N_A, N_B, N_C, N_A-1, N_B+1, N_C,N_A+1, N_B-1, N_C, 2000*dh) 
      else:
        dh = 9999
      dh_dict['%d_%d_%d'%(N_A, N_B, N_C)] = 2000 * dh
  return dh_dict

def get_Satoms(structure, atom_1,  index_list, n_neigh=4):
  from pylada.crystal import neighbors
  atom_indexes = []
  for id in index_list:
    atom = structure[id]
    p = atom.pos
    neighs = [n for n in neighbors(structure, n_neigh, p)]
    for n in neighs:
      if n[0].type == atom_1:
        index = list(structure).index(n[0])
        if index not in atom_indexes: 
          atom_indexes.append(index)
  return atom_indexes

def remove_atoms(structure, atoms_0= ['S',], atoms = ['Cu','Sn', 'Zn'], confs = [[3,1,0],],  remove_cluster=True, remove_metal=True):
    from pylada.crystal import neighbors
    structure_tmp = structure.copy()

    
    if remove_cluster:  
        atoms_notto_remove = []
        
        index_list = range(len(structure))
        for id in index_list:
            atom = structure[id]
            if atom.type in atoms_0:
                atoms_tmp = [atom]
                p = atom.pos
                neighs = [n for n in neighbors(structure, 4, p)]
                count = [0,0,0]
                for n in neighs:
                    atoms_tmp.append(n[0])
                    for i in range(3):
                      if n[0].type == atoms[i]:
                        count[i]  +=1
                    
                if count in confs:
                    atoms_notto_remove += atoms_tmp
                else:
                    pass
          
        N = len(structure)
        for i in range(N)[::-1]:
            atom = structure[i]
            if atom in atoms_notto_remove:
                pass
            else:
                structure_tmp.pop(i)
    if remove_metal:
        N = len(structure_tmp)
        for i in range(N)[::-1]:
            if structure_tmp[i].type not in atoms_0:
                structure_tmp.pop(i)
    return structure_tmp

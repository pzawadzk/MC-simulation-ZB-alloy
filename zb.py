
def  czts_stanite(name = ''):
  from pylada.crystal import Structure
  lattice = Structure( 2.477644960741404, 4.955156048010834, 0.000043820927176,
                   2.477648025452657, 0.000043821034886, 4.955159889217687,
                   4.985611280825058, 0.000045591817839, 0.000046761421437, 
                  scale  = 1.1000000000000001, name   = 'czts_kesterite_%s' %name)\
  .add_atom(2.477644711011844, 4.955207407631495, 2.492811889655014, 'Cu')\
  .add_atom(4.955202009303576, 2.477647055882982, 2.492810260184329, 'Cu')\
  .add_atom(2.477598363156763, 2.477597983355954, 0.000049167464818, 'Zn')\
  .add_atom(2.477646976241861, 2.477649975716003, 4.985602199906752, 'Sn')\
  .add_atom(3.684748302216236, 3.684764581470720, 3.665214351862760,  'S')\
  .add_atom(6.225737629203348, 6.225730692974708, 3.665244812907645,  'S')\
  .add_atom(3.748057631800364, 1.207133684391553, 1.320517482113959,  'S')\
  .add_atom(1.207143761314643, 3.748073706271233, 1.320518195542458,  'S')
# End of lattice definition.
  return  lattice

def  czts_kesterite(name = ''):
  from pylada.crystal import Structure
  lattice = Structure( 5.46435, -1.77e-07, 0,
           1.77e-07, 5.46435, 0,
           0, 0, 10.9105,
           scale=1,name   = 'czts_kesterite_%s' %name)\
  .add_atom(0, 0, 0, 'Cu')\
  .add_atom(2.73218, 2.73218, 5.45523, 'Cu')\
  .add_atom(-8.85e-08, 2.73218, 2.72761, 'Cu')\
  .add_atom(2.73218, 8.85e-08, 8.18284, 'Cu')\
  .add_atom(-8.85e-08, 2.73218, 8.18284, 'Zn')\
  .add_atom(2.73218, 8.85e-08, 2.72761, 'Zn')\
  .add_atom(0, 0, 5.45523, 'Sn')\
  .add_atom(2.73218, 2.73218, 0, 'Sn')\
  .add_atom(4.0613, 4.02738, 4.04931, 'S')\
  .add_atom(1.40305, 1.43697, 4.04931, 'S')\
  .add_atom(1.43697, 4.0613, 6.86114, 'S')\
  .add_atom(4.02738, 1.40305, 6.86114, 'S')\
  .add_atom(1.32913, 1.29521, 9.50454, 'S')\
  .add_atom(4.13522, 4.16914, 9.50454, 'S')\
  .add_atom(4.16914, 1.32913, 1.40592, 'S')\
  .add_atom(1.29521, 4.13522, 1.40592, 'S')
# End of lattice definition.
  return  lattice


def  zb(name = ''):
  from pylada.crystal import Structure
  lattice = Structure(0.00000, 2.50000, 2.50000,
                      2.50000, 0.00000, 2.50000,
                      2.50000, 2.50000, 0.00000, 
                    scale  = 1.0000000000000001, name   = 'zb_%s' %name)\
  .add_atom(0, 0, 0, 'Cu')\
  .add_atom(1.25, 1.25, 1.25, 'S')
# End of lattice definition.
  return  lattice


#
# Create a ASCII file of the following structure given the ouput of a PFLOTRAN
# simulation with ugi mesh:
#
# the number of BC1 connections
# the number of BC2 connections
# the number of internal connections
# list of the TRUE face area (id1 id2 area)
#
# TRUE face area is the product of projected face area and Orthogonality error

import h5py
import numpy as np

def check_connections_ids_and_area(f_in, struct=True, ugi=False):
  #open pflotran output
  f = h5py.File(f_in,'r')
  if struct:
    areas = np.array(f["Time:  5.00000E+01 y/Face_Area [m^2]"])
    connections = np.array(f["Connection Ids"])
  elif ugi:
    areas = np.array(f["   1 Time  5.00000E+01 y/Face Area [m^2]"])
    connections = np.array(f["Domain/Connection Ids"])
    cos_angle_non_orth = np.array(f["   1 Time  5.00000E+01 y/Face Non Orthogonality Angle"])
    areas = areas / (1 - cos_angle_non_orth)
  else: #uge
    areas = np.array(f["   1 Time  5.00000E+01 y/Face Area [m^2]"])
    connections = np.array(f["Connection Ids"])
  f.close()
  
  n_connections_internal = 0
  n_connections_bc1, n_connections_bc2 = 0, 0
  for x in connections:
    if (x[0] < 0 or x[1] < 0): 
      if (x[0] == -1 or x[1] == -1): n_connections_bc1 += 1
      if (x[0] == -2 or x[1] == -2): n_connections_bc2 += 1
    else: n_connections_internal += 1
  
  #write connections #
  if struct: out_name = "check_struct.txt"
  else: out_name = "check_unstruct.txt"
  out = open(out_name,'w')
  out.write(str(n_connections_bc1)+'\n')
  out.write(str(n_connections_bc2)+'\n')
  out.write(str(n_connections_internal)+'\n')
  
  #write connection detail
  #sort connections
  connections = np.sort(connections)
  ind = np.lexsort((connections[:,1],connections[:,0]))
  print(connections)
  print(connections[ind])
  for con,area in zip(connections[ind],areas[ind]):
    out.write("{} {} {:.6e}\n".format(con[0],con[1],area))

  return

  
if __name__ == "__main__":
  check_connections_ids_and_area("face_area_ugi.h5", struct=False, ugi=True)

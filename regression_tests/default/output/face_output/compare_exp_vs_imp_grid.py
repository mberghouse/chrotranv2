#
# Compare the output of PRINT_CONNECTION_ID and FACE_AREA outputed by PFLOTRAN
# and the connection and face area from the same grid in explicit format
#

import h5py
import numpy as np

def check_connections_ids_and_area_unstructured(verbose=True):
  #open pflotran output
  f = h5py.File("face_area.h5",'r')
  connections = np.array(f["Domain/Connection Ids"])
  areas = np.array(f["   1 Time  5.00000E+01 y/Face Area [m^2]"])
  cos_angle_non_orth = np.array(f["   1 Time  5.00000E+01 y/Face Non Orthogonality Angle"])
  f.close()
  #for i,a in enumerate(connections): print(a,cos_angle_non_orth[i])
  
  n_connections_internal = 0
  n_connections_bc1, n_connections_bc2 = 0, 0
  for x in connections:
    if (x[0] < 0 or x[1] < 0): 
      if (x[0] == -1 or x[1] == -1): n_connections_bc1 += 1
      if (x[0] == -2 or x[1] == -2): n_connections_bc2 += 1
    else: n_connections_internal += 1
  
  
  #open reference explicit grid
  src = open("ref.uge",'r')
  n_cells = int(src.readline().split()[-1])
  for i in range(n_cells): src.readline()
  n_connections_internal_ref = int(src.readline().split()[-1])
  area_dict = {}
  for i in range(n_connections_internal_ref):
    line = src.readline().split()
    con = tuple([int(x) for x in line[:2]])
    area = float(line[-1])
    area_dict[con] = area
  src.close()
  
  #bc1
  src = open("ref_BC1.ex",'r')
  n_connections_bc1_ref = int(src.readline().split()[-1])
  for i in range(n_connections_bc1_ref):
    line = src.readline().split()
    con = tuple([int(line[0]),-1])
    area = float(line[-1])
    area_dict[con] = area
  src.close()
  #bc2
  src = open("ref_BC2.ex",'r')
  n_connections_bc2_ref = int(src.readline().split()[-1])
  for i in range(n_connections_bc2_ref):
    line = src.readline().split()
    con = tuple([int(line[0]),-2])
    area = float(line[-1])
    area_dict[con] = area
  src.close()
  
  
  #compare
  error = False
  if (n_connections_internal != n_connections_internal_ref or
      n_connections_bc1 != n_connections_bc1_ref or
      n_connections_bc2 != n_connections_bc2_ref): error = True
      
  if verbose:
    print(f"\nNumber of internal connections: {n_connections_internal} / ",end="")
    print(f"{n_connections_internal_ref}")
    print(f"Number of BC1 connections: {n_connections_bc1} / {n_connections_bc1_ref}")
    print(f"Number of BC2 connections: {n_connections_bc2} / {n_connections_bc2_ref}")
  
  area_match = 0
  if verbose: print("Connection\tArea\tRef\tError")
  for con,area in area_dict.items():
    found = False
    if verbose: print(f"{con}\t",end="")
    i = 0
    for x in connections:
      if (x[0] == con[0] and x[1] == con[1]) or \
         (x[1] == con[0] and x[0] == con[1]):
        found = True
        if not (x[0] < 0 or x[1] < 0): 
          true_area = areas[i] / (1-cos_angle_non_orth[i])
        else: true_area = areas[i]
        err_rel_area = abs(true_area-area) / area
        if verbose:
          print(f"{true_area:.6f}\t{area:.6f}\t{err_rel_area:6f}",end="")
          if err_rel_area < 1e-6: print("  OK")
          else: print("  X")
        if err_rel_area > 1e-6: error = True  
        break
      else: i += 1
    if not found:
      error = True
      if verbose: print(f"Connection not found!")
  
  if verbose:
    if error: print("\nError found !!!")
    else: print("\nTest passed")
  return error



def check_connections_ids_and_area_structured(verbose=True):
  #open pflotran output
  f = h5py.File("face_area_struct.h5",'r')
  connections = np.array(f["Connection Ids"])
  areas = np.array(f["Time:  5.00000E+01 y/Face_Area [m^2]"])
  #cos_angle_non_orth = np.array(f["   1 Time  5.00000E+01 y/Face Non Orthogonality Angle"])
  f.close()
  #for i,a in enumerate(connections): print(a,cos_angle_non_orth[i])
  
  n_connections_internal = 0
  n_connections_bc1, n_connections_bc2 = 0, 0
  for x in connections:
    if (x[0] < 0 or x[1] < 0): 
      if (x[0] == -1 or x[1] == -1): n_connections_bc1 += 1
      if (x[0] == -2 or x[1] == -2): n_connections_bc2 += 1
    else: n_connections_internal += 1
  
  
  #open reference explicit grid
  src = open("ref_struct.uge",'r')
  n_cells = int(src.readline().split()[-1])
  for i in range(n_cells): src.readline()
  n_connections_total = int(src.readline().split()[-1])
  area_dict = {}
  n_connections_bc1_ref = 0
  n_connections_bc2_ref = 0
  n_connections_internal_ref = 0
  for i in range(n_connections_total):
    line = src.readline().split()
    con = tuple([int(x) for x in line[:2]])
    area = float(line[-1])
    if (con[1] == -1): n_connections_bc1_ref += 1
    elif (con[1] == -2): n_connections_bc2_ref += 1
    else: n_connections_internal_ref += 1
    area_dict[con] = area
  src.close()
  
  
  #compare
  error = False
  if (n_connections_internal != n_connections_internal_ref or
      n_connections_bc1 != n_connections_bc1_ref or
      n_connections_bc2 != n_connections_bc2_ref): error = True
      
  if verbose:
    print(f"\nNumber of internal connections: {n_connections_internal} / ",end="")
    print(f"{n_connections_internal_ref}")
    print(f"Number of BC1 connections: {n_connections_bc1} / {n_connections_bc1_ref}")
    print(f"Number of BC2 connections: {n_connections_bc2} / {n_connections_bc2_ref}")
  
  area_match = 0
  if verbose: print("Connection\tArea\tRef\tError")
  for con,area in area_dict.items():
    found = False
    if verbose: print(f"{con}\t",end="")
    i = 0
    for x in connections:
      if (x[0] == con[0] and x[1] == con[1]) or \
         (x[1] == con[0] and x[0] == con[1]):
        found = True
        true_area = areas[i]
        err_rel_area = abs(true_area-area) / area
        if verbose:
          print(f"{true_area:.6f}\t{area:.6f}\t{err_rel_area:6f}",end="")
          if err_rel_area < 1e-6: print("  OK")
          else: print("  X")
        if err_rel_area > 1e-6: error = True  
        break
      else: i += 1
    if not found:
      error = True
      if verbose: print(f"Connection not found!")
  
  if verbose:
    if error: print("\nError found !!!")
    else: print("\nTest passed")
  return error
  

  
if __name__ == "__main__":
  err = check_connections_ids_and_area_unstructured(verbose = True)
  err2 = check_connections_ids_and_area_structured(verbose = True)
  if err or err2: 
    print("!!! Some test errored !!!")
    exit(-1)
  else: 
    print("All tests passed")
    exit(0)

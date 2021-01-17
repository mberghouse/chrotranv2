import numpy as np
import h5py

def compare_jacobian(f,ref="sensitivity_pressure.ref"):
  #reference file
  J_ref = np.genfromtxt(ref, skip_header=8, skip_footer=2)
  rows_ref, cols_ref, datas_ref = J_ref[:,0], J_ref[:,1], J_ref[:,2]
  
  #get hdf5 sensitivity
  src = h5py.File(f,'r')
  rows = np.array(src["Mat Structure/Row Indices"])
  cols = np.array(src["Mat Structure/Column Indices"])
  datas = np.array(src['Time:  1.00000E+02 y/Pressure []'])
  
  if len(rows_ref) != len(rows): 
    print("ERROR: Matrices does not have the same structure")
  
  error = False
  print("I\tJ\tValue")
  for i in range(len(rows_ref)):
    found = False
    row_ref = rows_ref[i]
    col_ref = cols_ref[i]
    print(f"{int(row_ref)}\t{int(col_ref)}\t",end="")
    for j in range(len(rows)):
      x = rows[j]
      y = cols[j]
      if row_ref == x and col_ref == y:
        found = True
        print(f"{datas[j]:.6E} vs {datas_ref[i]:.6E}\t",end='')
        if abs(1-datas[j]/datas_ref[i]) < 1e-6:
          print("OK")
        else:
          error = True
          print("X")
        break
    if not found:
      error = True
      print("Element not found")
  
  if error: 
    print("\nSome error occured")
    return 1
  else: 
    print("\nMatrices are the same")
    return 0

if __name__ == "__main__":
  err = compare_jacobian("imp_128-sensitivity-flow.h5")
  exit(err)

  

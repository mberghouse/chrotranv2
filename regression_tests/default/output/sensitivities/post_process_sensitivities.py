import numpy as np
import h5py


def post_process(f):

  src = h5py.File(f+"-sensitivity-flow.h5",'r')
  
  datasets = ["Permeability","Pressure"] 
  rows = np.array(src["Mat Structure/Row Indices"])
  cols = np.array(src["Mat Structure/Column Indices"])
  indexes = np.zeros((len(rows),2),dtype='i8')
  count = 0
  for i,j in zip(rows,cols):
    indexes[count] = i,j
    count += 1
  order = np.lexsort((indexes[:,1],indexes[:,0]))

  for name in datasets:
    data = np.array(src[f'Time:  1.00000E+02 y/{name} []'])
    out = open(f"sensitivity-{name.lower()}.csv",'w')
    for i in order:
      out.write(f"{rows[i]} {indexes[i][1]}  {data[i]:.6e}\n")
    out.close()
  
  src.close()
  
  

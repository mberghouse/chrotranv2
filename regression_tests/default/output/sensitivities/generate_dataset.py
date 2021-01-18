import numpy as np
import h5py

if __name__ == "__main__":
  #read permeability from csv
  perm_data = np.genfromtxt("permeability_field.txt", comments='#')
  cell_ids, perm_field = perm_data[:,0], perm_data[:,1]
  #create dataset
  out = h5py.File("permeability.h5", 'w')
  out.create_dataset("Cell Ids", data=cell_ids)
  out.create_dataset("permeability", data=perm_field)
  out.close()
    

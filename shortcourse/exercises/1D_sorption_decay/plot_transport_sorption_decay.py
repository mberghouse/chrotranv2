import sys
import os
try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  try:
    pflotran_dir = '../../'
  except KeyError:
    print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
    sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from plot_aux import *

f = plt.figure(figsize=(16,6))
f.suptitle("Solute Transport, Sorption and Decay",fontsize=16)

plt.subplot(1,2,1)
plt.title('Concentration Breakthrough @ 49.5 Meters')
plt.xlabel('Concentration [M]')
plt.ylabel('Time [y]')
filenames = ['transport_sorption_decay-obs-0.pft']
columns = np.arange(2,8)
minval,maxval = plot_results(plt,filenames,columns)


plt.subplot(1,2,2)
plt.title('Concentration Profile @ 12.5 Years')
plt.xlabel('Concentration [M]')
plt.ylabel('X [m]')
filenames = ['transport_sorption_decay-002.tec']
columns = np.arange(4,10)
minval,maxval = plot_results(plt,filenames,columns)

f.subplots_adjust(hspace=0.2,wspace=0.40,
                  bottom=.12,top=.85,
                  left=.08,right=.92)

plt.show()

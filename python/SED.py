import numpy as np
import matplotlib.pyplot as plt
import os

path = '/scratch/ntf229/880Msun_split/SKIRT/'

plt.figure(figsize=(10,8))

wave = np.load(path+'wave.npy') # units of microns
wave = wave * 1e4 # convert to Angstroms
spec = np.load(path+'spec.npy') # units of Jy
spec = spec / 3631 # convert to maggies

mask = (wave >= 1e3) & (wave <= 1e7)

plt.plot(wave[mask], spec[mask], color='blue',alpha=1)
plt.xlabel('Wavelength ('+r'$\AA$'+')', fontsize=16)
plt.ylabel('Flux (Maggies)', fontsize=16)

plt.xscale('log')
plt.yscale('log')

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(path+'SED.png', dpi=300)
plt.close()


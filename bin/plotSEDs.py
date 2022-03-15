import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--index") # index of snapshot (see snaps below)                                                   $
args = parser.parse_args()
index = int(args.index)

snaps = ["snapdir_570","snapdir_576","snapdir_581","snapdir_586","snapdir_591","snapdir_596","snapdir_601","snapdir_606","snapdir_611","snapdir_616","snapdir_621","snapdir_626","snapdir_631","snapdir_636","snapdir_641","snapdir_646","snapdir_651","snapdir_656","snapdir_661","snapdir_666","snapdir_671","snapdir_676","snapdir_681","snapdir_686","snapdir_691","snapdir_696"]

SKIRTpath = '/scratch/ntf229/880Msun_split/SKIRT/'+snaps[index]+'/'
plotPath = '/scratch/ntf229/880Msun_split/SEDs/'

os.system('mkdir -p '+plotPath)

plt.figure(figsize=(10,8))

wave = np.load(SKIRTpath+'wave.npy') # units of microns
wave = wave * 1e4 # convert to Angstroms
spec = np.load(SKIRTpath+'spec.npy') # units of Jy
spec = spec / 3631 # convert to maggies

mask = (wave >= 1e3) & (wave <= 1e7)

plt.plot(wave[mask], spec[mask], color='blue',alpha=1)
plt.xlabel('Wavelength ('+r'$\AA$'+')', fontsize=16)
plt.ylabel('Flux (Maggies)', fontsize=16)

plt.xscale('log')
plt.yscale('log')

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(plotPath+snaps[index]+'_SED.png', dpi=300)
plt.close()


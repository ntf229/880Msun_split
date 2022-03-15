import numpy as np
import matplotlib.pyplot as plt
import os

snaps = ["snapdir_570","snapdir_576","snapdir_581","snapdir_586","snapdir_591","snapdir_596","snapdir_601","snapdir_606","snapdir_611","snapdir_616","snapdir_621","snapdir_626","snapdir_631","snapdir_636","snapdir_641","snapdir_646","snapdir_651","snapdir_656","snapdir_661","snapdir_666","snapdir_671","snapdir_676","snapdir_681","snapdir_686","snapdir_691","snapdir_696"]
plotPath = '/scratch/ntf229/880Msun_split/allSEDs/'
os.system('mkdir -p '+plotPath)
plt.figure(figsize=(10,8))
colors = plt.cm.rainbow(np.linspace(0, 1, len(snaps)-1))

for i in range(len(snaps)):
    SKIRTpath = '/scratch/ntf229/880Msun_split/SKIRT/'+snaps[i]+'/'
    wave = np.load(SKIRTpath+'wave.npy') # units of microns
    wave = wave * 1e4 # convert to Angstroms
    spec = np.load(SKIRTpath+'spec.npy') # units of Jy
    spec = spec / 3631 # convert to maggies
    mask = (wave >= 1e3) & (wave <= 1e7)
    if i == 0:
        plt.plot(wave[mask], spec[mask], color='k', alpha=1, label='parent') # parent
    else:
        plt.plot(wave[mask], spec[mask], color=colors[i-1], alpha=1)

plt.xlabel('Wavelength ('+r'$\AA$'+')', fontsize=16)
plt.ylabel('Flux (Maggies)', fontsize=16)

plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=16)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(plotPath+'allSEDs.png', dpi=300)
plt.close()


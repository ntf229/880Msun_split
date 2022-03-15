import numpy as np
import matplotlib.pyplot as plt
import os

snaps = ["snapdir_570","snapdir_576","snapdir_581","snapdir_586","snapdir_591","snapdir_596","snapdir_601","snapdir_606","snapdir_611","snapdir_616","snapdir_621","snapdir_626","snapdir_631","snapdir_636","snapdir_641","snapdir_646","snapdir_651","snapdir_656","snapdir_661","snapdir_666","snapdir_671","snapdir_676","snapdir_681","snapdir_686","snapdir_691","snapdir_696"]
plotPath = '/scratch/ntf229/880Msun_split/Totals/'
os.system('mkdir -p '+plotPath)
plt.figure(figsize=(10,8))
colors = plt.cm.rainbow(np.linspace(0, 1, len(snaps)-1))

for i in range(len(snaps)):
    textPath = '/scratch/ntf229/880Msun_split/TextFiles/'+snaps[i]+'/'
    stars = np.loadtxt(textPath+'stars.txt')
    starMasses = stars[:,7] # units of M_sun
    totalStellarMass = np.sum(starMasses)
    starAges = stars[:,9] # units of years
    ageMask = starAges < 1e8 # younger than 100 Myrs
    youngStars = np.loadtxt(textPath+'youngStars.txt')
    if len(youngStars) == 0:
        youngTotalMass = 0.
    else:
        youngMass = youngStars[:,7] * 1e7
        youngTotalMass = np.sum(youngMass)
    totalStellarMass += youngTotalMass
    SFR = (np.sum(starMasses[ageMask]) + youngTotalMass) / 1e8 # solar masses per year
    if i == 0:
        parentMass = totalStellarMass
        parentSFR = SFR
        print('saving parent values')
        print('parent mass:',parentMass)
        print('parent SFR:',parentSFR)
        #plt.scatter(totalStellarMass, SFR, color='k', alpha=1, label='parent', s=15 # parent
    else:
        plt.scatter(totalStellarMass, SFR, color=colors[int(i-1)], alpha=1, s=10)

plt.scatter(parentMass, parentSFR, color='k', alpha=1, label='parent', s=15) # parent
plt.xlabel('Total Stellar Mass ('+r'$M_{\odot}$'+')', fontsize=16)
plt.ylabel('SFR ('+r'$M_{\odot}$'+' / year)', fontsize=16)

plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=16)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(plotPath+'totals.png', dpi=300)
plt.close()


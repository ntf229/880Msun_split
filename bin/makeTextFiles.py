# make text files for SKIRT input

import numpy as np
import os
import matplotlib.pyplot as plt
import h5py

filePath = '/scratch/ntf229/880Msun_split/uv_background/m12i_res7100_uvb-late/1Myr_res880/output/snapdir_696/'
textPath = '/scratch/ntf229/880Msun_split/TextFiles/'
plotPath = '/home/ntf229/CCA/880Msun_split/plots/'

# Create directories if they don't already exist
os.system('mkdir -p '+textPath)
os.system('mkdir -p '+plotPath)

# Spatial cuts (in pc)
xAbsMin = 4.180e7 # x,y,z length: 0.005e7 pc = 5e4 pc = 50 kpc 
xAbsMax = 4.185e7
yAbsMin = 4.4125e7
yAbsMax = 4.4175e7
zAbsMin = 4.625e7
zAbsMax = 4.630e7

# constants
proton_mass = 8.4089382e-58 # in solar masses
k_Boltzmann = 6.94169e-60 # in km**2 M_sun s**-2 K**-1
gamma = 5./3.
XH = 0.76 
y_helium = (1.0-XH)/(4.0*XH)  # from xiangchang ma physics.py

f = h5py.File(filePath+'snapshot_696.0.hdf5', 'r')
numFiles = f['Header'].attrs['NumFilesPerSnapshot']
hubble = f['Header'].attrs['HubbleParam']
hinv = 1.0 / hubble
OmegaM = 0.272 

for i in range(numFiles):
        
    # load data from current file
    f = h5py.File(filePath+'snapshot_696.'+str(i)+'.hdf5', 'r')
    stars = f['PartType4']
    gas = f['PartType0']
    
    # stars
    x_pos = stars['Coordinates'][:,0] * hinv * 1e3 # in pc
    y_pos = stars['Coordinates'][:,1] * hinv * 1e3
    z_pos = stars['Coordinates'][:,2] * hinv * 1e3
    x_vel = stars['Velocities'][:,0] # in km/s
    y_vel = stars['Velocities'][:,1]
    z_vel = stars['Velocities'][:,2]
    mass = stars['Masses'] * hinv * 1e10 # in solar masses
    smooth = np.zeros(len(mass)) + 4 # assume all stars have 4pc smoothing length
    metals = stars['Metallicity'][:,0] # in solar units (mass fraction)
    z = 1./stars['StellarFormationTime'][:] - 1
    x = OmegaM / (1.0-OmegaM) * (1+z)**3
    t = (2.0/(3.0*np.sqrt(1.0-OmegaM))) * np.log(np.sqrt(x)/(-1.0+np.sqrt(1.+x)))
    t *= (13.777*(0.71/hubble)) # in Gyr
    t *= 1e9 # convert to years
    age = 1.38e10 - t
    
    # gas
    x_pos_gas = gas['Coordinates'][:,0] * hinv * 1e3 # in pc
    y_pos_gas = gas['Coordinates'][:,1] * hinv * 1e3
    z_pos_gas = gas['Coordinates'][:,2] * hinv * 1e3
    smooth_gas = gas['SmoothingLength'] * hinv * 1e3 # in pc
    mass_gas = gas['Masses'] * hinv * 1e10 # in solar masses
    metals_gas = gas['Metallicity'][:,0] # in solar units (mass fraction)
    ElectronAbundance = np.asarray(gas['ElectronAbundance'][:])
    helium_mass_fraction = gas['Metallicity'][:,1]
    InternalEnergy = gas['InternalEnergy'][:] #  * 1e10
    mu = np.asarray((1.+4.*y_helium)/(1.+y_helium+ElectronAbundance))
    mean_molecular_weight = np.zeros(len(mu))
    temp_gas = np.zeros(len(mu))
    for j in range(len(mu)):
        mean_molecular_weight[j] = mu[j]*proton_mass
        temp_gas[j] = mean_molecular_weight[j] * (gamma-1) * InternalEnergy[j] / k_Boltzmann

    # append data
    if i == 0:
        full_x_pos = x_pos 
        full_y_pos = y_pos
        full_z_pos = z_pos
        full_smooth = smooth
        full_x_vel = x_vel
        full_y_vel = y_vel
        full_z_vel = z_vel
        full_mass = mass
        full_metals = metals
        full_age = age
        full_x_pos_gas = x_pos_gas
        full_y_pos_gas = y_pos_gas
        full_z_pos_gas = z_pos_gas
        full_smooth_gas = smooth_gas
        full_mass_gas = mass_gas
        full_metals_gas = metals_gas
        full_temp_gas = temp_gas
    else:
        full_x_pos = np.append(full_x_pos, x_pos)
        full_y_pos = np.append(full_y_pos, y_pos)
        full_z_pos = np.append(full_z_pos, z_pos)
        full_smooth = np.append(full_smooth, smooth)
        full_x_vel = np.append(full_x_vel, x_vel)
        full_y_vel = np.append(full_y_vel, y_vel)
        full_z_vel = np.append(full_z_vel, z_vel)
        full_mass = np.append(full_mass, mass)
        full_metals = np.append(full_metals, metals)
        full_age = np.append(full_age, age)
        full_x_pos_gas = np.append(full_x_pos_gas, x_pos_gas)
        full_y_pos_gas = np.append(full_y_pos_gas, y_pos_gas)
        full_z_pos_gas = np.append(full_z_pos_gas, z_pos_gas)
        full_smooth_gas = np.append(full_smooth_gas, smooth_gas)
        full_mass_gas = np.append(full_mass_gas, mass_gas)
        full_metals_gas = np.append(full_metals_gas, metals_gas)
        full_temp_gas = np.append(full_temp_gas, temp_gas)

# Make cuts in the data based on spatial position
starIndex = [] # indices of stars to be removed
gasIndex = [] # indices of gas to be removed

for i in range(len(full_x_pos)):
    if (full_x_pos[i] < xAbsMin) or (full_x_pos[i] > xAbsMax):
        starIndex.append(i)
    elif (full_y_pos[i] < yAbsMin) or (full_y_pos[i] > yAbsMax):
        starIndex.append(i)
    elif (full_z_pos[i] < zAbsMin) or (full_z_pos[i] > zAbsMax):
        starIndex.append(i)

for i in range(len(full_x_pos_gas)):
    if (full_x_pos_gas[i] < xAbsMin) or (full_x_pos_gas[i] > xAbsMax):
        gasIndex.append(i)
    elif (full_y_pos_gas[i] < yAbsMin) or (full_y_pos_gas[i] > yAbsMax):
        gasIndex.append(i)
    elif (full_z_pos_gas[i] < zAbsMin) or (full_z_pos_gas[i] > zAbsMax):
        gasIndex.append(i)

cut_x_pos = np.float32(np.delete(full_x_pos,starIndex))
cut_y_pos = np.float32(np.delete(full_y_pos,starIndex))
cut_z_pos = np.float32(np.delete(full_z_pos,starIndex))
cut_mass = np.float32(np.delete(full_mass,starIndex))
cut_metals = np.float32(np.delete(full_metals,starIndex) )
cut_age = np.float32(np.delete(full_age,starIndex))
cut_x_vel = np.float32(np.delete(full_x_vel,starIndex))
cut_y_vel = np.float32(np.delete(full_y_vel,starIndex))
cut_z_vel = np.float32(np.delete(full_z_vel,starIndex))
cut_smooth = np.float32(np.delete(full_smooth,starIndex))

cut_x_pos_gas = np.float32(np.delete(full_x_pos_gas,gasIndex))
cut_y_pos_gas = np.float32(np.delete(full_y_pos_gas,gasIndex))
cut_z_pos_gas = np.float32(np.delete(full_z_pos_gas,gasIndex))
cut_smooth_gas = np.float32(np.delete(full_smooth_gas,gasIndex))
cut_mass_gas = np.float32(np.delete(full_mass_gas,gasIndex))
cut_metals_gas = np.float32(np.delete(full_metals_gas,gasIndex))
cut_temp_gas = np.float32(np.delete(full_temp_gas,gasIndex))

# Shift all particles to be centered at x,y,z = 0
xCenter = (xAbsMax + xAbsMin)/2
yCenter = (yAbsMax + yAbsMin)/2
zCenter = (zAbsMax + zAbsMin)/2

cut_x_pos = cut_x_pos - xCenter
cut_y_pos = cut_y_pos - yCenter
cut_z_pos = cut_z_pos - zCenter

cut_x_pos_gas = cut_x_pos_gas - xCenter
cut_y_pos_gas = cut_y_pos_gas - yCenter
cut_z_pos_gas = cut_z_pos_gas - zCenter 

print('Number of star particles:',len(cut_x_pos))
print('Number of gas particles:',len(cut_x_pos_gas))

star_header = 'Column 1: position x (pc)\nColumn 2: position y (pc)\nColumn 3: position z (pc)\nColumn 4: smoothing length (pc)\nColumn 5: v_x (km/s)\nColumn 6: v_y (km/s)\nColumn 7: v_z (km/s)\nColumn 8: mass (Msun)\nColumn 9: metallicity ()\nColumn 10: age (yr)'
gas_header = 'Column 1: position x (pc)\nColumn 2: position y (pc)\nColumn 3: position z (pc)\nColumn 4: smoothing length (pc)\nColumn 5: mass (Msun)\nColumn 6: metallicity ()\nColumn 7: temperature (K)' 

np.savetxt(textPath+'stars.txt',np.float32(np.c_[cut_x_pos, cut_y_pos, cut_z_pos, cut_smooth, cut_x_vel, cut_y_vel, cut_z_vel, cut_mass, cut_metals, cut_age]),header=star_header)
np.savetxt(textPath+'gas.txt',np.float32(np.c_[cut_x_pos_gas, cut_y_pos_gas, cut_z_pos_gas, cut_smooth_gas, cut_mass_gas, cut_metals_gas, cut_temp_gas]),header=gas_header)


exit()

# Make plots

# SFH 
plt.figure(figsize=(10,8))
counts, bins = np.histogram(cut_age,bins=100,weights=cut_mass,density=False)
plt.hist(bins[:-1], bins, weights=counts)
plt.xlabel('Age (years)', fontsize=16)
plt.ylabel('Stellar Mass ('+r'$M_{\odot}$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.yscale('log')
plt.savefig(plotPath+'SFH.png',dpi=300)
plt.close()

# Gas Metallicity 
plt.figure(figsize=(10,8))
counts, bins = np.histogram(cut_metals_gas,bins=500,weights=cut_mass_gas,density=False)
plt.hist(bins[:-1], bins, weights=counts)
plt.xlabel('Metallicity', fontsize=16)
plt.ylabel('Gas Mass ('+r'$M_{\odot}$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.yscale('log')
plt.savefig(plotPath+'gasMetals.png',dpi=300)
plt.close()

# Gas Temperature Distribution
plt.figure(figsize=(10,8))
counts, bins = np.histogram(np.log10(cut_temp_gas),bins=500,weights=cut_mass_gas,density=False)
plt.hist(bins[:-1], bins, weights=counts)
plt.xlabel('Log(T)', fontsize=16)
plt.ylabel('Gas Mass ('+r'$M_{\odot}$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.yscale('log')
plt.savefig(plotPath+'gasTemp.png',dpi=300)
plt.close()

# gas x position
plt.figure(figsize=(10,8))
counts, bins = np.histogram(cut_x_pos_gas,bins=500,weights=cut_mass_gas,density=False)
plt.hist(bins[:-1], bins, weights=counts)
plt.xlabel('x (pc)', fontsize=16)
plt.ylabel('Gas Mass ('+r'$M_{\odot}$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.yscale('log')
plt.savefig(plotPath+'gasX.png',dpi=300)
plt.close()

# gas y position
plt.figure(figsize=(10,8))
counts, bins = np.histogram(cut_y_pos_gas,bins=500,weights=cut_mass_gas,density=False)
plt.hist(bins[:-1], bins, weights=counts)
plt.xlabel('y (pc)', fontsize=16)
plt.ylabel('Gas Mass ('+r'$M_{\odot}$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.yscale('log')
plt.savefig(plotPath+'gasY.png',dpi=300)
plt.close()

# gas z position
plt.figure(figsize=(10,8))
counts, bins = np.histogram(cut_z_pos_gas,bins=500,weights=cut_mass_gas,density=False)
plt.hist(bins[:-1], bins, weights=counts)
plt.xlabel('z (pc)', fontsize=16)
plt.ylabel('Gas Mass ('+r'$M_{\odot}$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.yscale('log')
plt.savefig(plotPath+'gasZ.png',dpi=300)
plt.close()

# stars x position
plt.figure(figsize=(10,8))
counts, bins = np.histogram(cut_x_pos,bins=500,weights=cut_mass,density=False)
plt.hist(bins[:-1], bins, weights=counts)
plt.xlabel('x (pc)', fontsize=16)
plt.ylabel('Stellar Mass ('+r'$M_{\odot}$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.yscale('log')
plt.savefig(plotPath+'starsX.png',dpi=300)
plt.close()

# stars y position
plt.figure(figsize=(10,8))
counts, bins = np.histogram(cut_y_pos,bins=500,weights=cut_mass,density=False)
plt.hist(bins[:-1], bins, weights=counts)
plt.xlabel('y (pc)', fontsize=16)
plt.ylabel('Stellar Mass ('+r'$M_{\odot}$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.yscale('log')
plt.savefig(plotPath+'starsY.png',dpi=300)
plt.close()

# stars z position
plt.figure(figsize=(10,8))
counts, bins = np.histogram(cut_z_pos,bins=500,weights=cut_mass,density=False)
plt.hist(bins[:-1], bins, weights=counts)
plt.xlabel('z (pc)', fontsize=16)
plt.ylabel('Stellar Mass ('+r'$M_{\odot}$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.yscale('log')
plt.savefig(plotPath+'starsZ.png',dpi=300)
plt.close()






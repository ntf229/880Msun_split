# make text files for SKIRT input

import numpy as np
import os
import matplotlib.pyplot as plt
import h5py

filePath = '/scratch/ntf229/880Msun_split/uv_background/m12i_res7100_uvb-late/1Myr_res880/output/snapdir_696/'

storePath = '/scratch/ntf229/880Msun_split/resources/particles/'

# Create directories if they don't already exist
os.system('mkdir -p '+storePath)

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


np.save(storePath+'starIndex.npy', starIndex)
np.save(storePath+'gasIndex.npy', gasIndex)

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

center = [np.average(cut_x_pos, weights=cut_mass), np.average(cut_y_pos, weights=cut_mass), np.average(cut_z_pos, weights=cut_mass)] 
np.save(storePath+'center.npy', center)

print('done')

exit()


cut_x_pos_gas = np.float32(np.delete(full_x_pos_gas,gasIndex))
cut_y_pos_gas = np.float32(np.delete(full_y_pos_gas,gasIndex))
cut_z_pos_gas = np.float32(np.delete(full_z_pos_gas,gasIndex))
cut_smooth_gas = np.float32(np.delete(full_smooth_gas,gasIndex))
cut_mass_gas = np.float32(np.delete(full_mass_gas,gasIndex))
cut_metals_gas = np.float32(np.delete(full_metals_gas,gasIndex))
cut_temp_gas = np.float32(np.delete(full_temp_gas,gasIndex))


# make text files for SKIRT input

import numpy as np
import os
import matplotlib.pyplot as plt
import h5py
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--index") # index of snapshot (see snaps below)                                                                                        $
args = parser.parse_args()
index = int(args.index)

snaps = ["snapdir_570","snapdir_576","snapdir_581","snapdir_586","snapdir_591","snapdir_596","snapdir_601","snapdir_606","snapdir_611","snapdir_616","snapdir_621","snapdir_626","snapdir_631","snapdir_636","snapdir_641","snapdir_646","snapdir_651","snapdir_656","snapdir_661","snapdir_666","snapdir_671","snapdir_676","snapdir_681","snapdir_686","snapdir_691","snapdir_696"]
snap_num = snaps[index].split('_')[1]

if index == 0: # parent snapshot
    filePath = '/scratch/ntf229/880Msun_split/uv_background/m12i_res7100_uvb-late/output/'+snaps[index]+'/'
else:
    filePath = '/scratch/ntf229/880Msun_split/uv_background/m12i_res7100_uvb-late/1Myr_res880/output/'+snaps[index]+'/'

textPath = '/scratch/ntf229/880Msun_split/TextFiles/'+snaps[index]+'/'

# Create directories if they don't already exist
os.system('mkdir -p '+textPath)

# Spatial cuts for snapshot 696 (in pc)
xAbsMin = 4.180e7 # x,y,z length: 0.005e7 pc = 5e4 pc = 50 kpc 
xAbsMax = 4.185e7
yAbsMin = 4.4125e7
yAbsMax = 4.4175e7
zAbsMin = 4.625e7
zAbsMax = 4.630e7

vel696 = np.asarray([-50.949375, 73.33526, 97.57783]) # average stellar velocity from snapshot 696 in km/s
vel696 *= 3.241e-14 # conver to pc/s

timeInt = 696. - float(snap_num) # time interval between current snapshot and snapshot 696
timeInt *= 3.154e13 # convert from Myrs to seconds

x_offset = -1. * vel696[0] * timeInt
y_offset = -1. * vel696[1] * timeInt
z_offset = -1. * vel696[2] * timeInt

xAbsMin += x_offset
xAbsMax += x_offset
yAbsMin += y_offset
yAbsMax += y_offset
zAbsMin += z_offset
zAbsMax += z_offset

xCenter = (xAbsMin + xAbsMax) / 2
yCenter = (yAbsMin + yAbsMax) / 2
zCenter = (zAbsMin + zAbsMax) / 2

# constants
proton_mass = 8.4089382e-58 # in solar masses
k_Boltzmann = 6.94169e-60 # in km**2 M_sun s**-2 K**-1
gamma = 5./3.
XH = 0.76 
y_helium = (1.0-XH)/(4.0*XH)  # from xiangchang ma physics.py

#f = h5py.File(filePath+'snapshot_696.0.hdf5', 'r')
f = h5py.File(filePath+'snapshot_'+snap_num+'.0.hdf5', 'r')
numFiles = f['Header'].attrs['NumFilesPerSnapshot']
hubble = f['Header'].attrs['HubbleParam']
hinv = 1.0 / hubble
OmegaM = 0.272 

for i in range(numFiles):
    
    # load data from current file
    #f = h5py.File(filePath+'snapshot_696.'+str(i)+'.hdf5', 'r')
    f = h5py.File(filePath+'snapshot_'+snap_num+'.'+str(i)+'.hdf5', 'r')
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

# make spatial cuts  
boxSize = 5.e4 # box side length in pc

xMin = full_x_pos > xAbsMin
xMax = full_x_pos[xMin] < xAbsMax
yMin = full_y_pos[xMin][xMax] > yAbsMin
yMax = full_y_pos[xMin][xMax][yMin] < yAbsMax
zMin = full_z_pos[xMin][xMax][yMin][yMax] > zAbsMin
zMax = full_z_pos[xMin][xMax][yMin][yMax][zMin] < zAbsMax

cut_x_pos = full_x_pos[xMin][xMax][yMin][yMax][zMin][zMax]
cut_y_pos = full_y_pos[xMin][xMax][yMin][yMax][zMin][zMax]
cut_z_pos = full_z_pos[xMin][xMax][yMin][yMax][zMin][zMax]
cut_smooth = full_smooth[xMin][xMax][yMin][yMax][zMin][zMax]
cut_x_vel = full_x_vel[xMin][xMax][yMin][yMax][zMin][zMax]
cut_y_vel = full_y_vel[xMin][xMax][yMin][yMax][zMin][zMax]
cut_z_vel = full_z_vel[xMin][xMax][yMin][yMax][zMin][zMax]
cut_mass = full_mass[xMin][xMax][yMin][yMax][zMin][zMax]
cut_metals = full_metals[xMin][xMax][yMin][yMax][zMin][zMax]
cut_age = full_age[xMin][xMax][yMin][yMax][zMin][zMax]

xMinGas = full_x_pos_gas > xAbsMin
xMaxGas = full_x_pos_gas[xMinGas] < xAbsMax
yMinGas = full_y_pos_gas[xMinGas][xMaxGas] > yAbsMin
yMaxGas = full_y_pos_gas[xMinGas][xMaxGas][yMinGas] < yAbsMax
zMinGas = full_z_pos_gas[xMinGas][xMaxGas][yMinGas][yMaxGas] > zAbsMin
zMaxGas = full_z_pos_gas[xMinGas][xMaxGas][yMinGas][yMaxGas][zMinGas] < zAbsMax

cut_x_pos_gas = full_x_pos_gas[xMinGas][xMaxGas][yMinGas][yMaxGas][zMinGas][zMaxGas]
cut_y_pos_gas = full_y_pos_gas[xMinGas][xMaxGas][yMinGas][yMaxGas][zMinGas][zMaxGas]
cut_z_pos_gas = full_z_pos_gas[xMinGas][xMaxGas][yMinGas][yMaxGas][zMinGas][zMaxGas]
cut_smooth_gas = full_smooth_gas[xMinGas][xMaxGas][yMinGas][yMaxGas][zMinGas][zMaxGas]
cut_mass_gas = full_mass_gas[xMinGas][xMaxGas][yMinGas][yMaxGas][zMinGas][zMaxGas]
cut_metals_gas = full_metals_gas[xMinGas][xMaxGas][yMinGas][yMaxGas][zMinGas][zMaxGas]
cut_temp_gas = full_temp_gas[xMinGas][xMaxGas][yMinGas][yMaxGas][zMinGas][zMaxGas]

# shift cut particles to center
cut_x_pos -= xCenter
cut_y_pos -= yCenter
cut_z_pos -= zCenter
cut_x_pos_gas -= xCenter
cut_y_pos_gas -= yCenter
cut_z_pos_gas -= zCenter

print('number of star particles:', len(cut_x_pos))
print('number of gas particles:', len(cut_x_pos_gas))

# create MAPPINGS-III particles
tau_clear = 2.e6 # assume tau_clear years clearing time (e-folding time)
k = 1.3806485*10**(-19) # boltzmann constant in cm**2 kg s**-2 K**-1
youngStarMask = cut_age < 1.e7 # Mask for young stars 
youngStarIndex = [] # Indices of young stars 
for i in range(len(cut_age)):
    if youngStarMask[i]:
        youngStarIndex.append(i)
print('number of young star particles:',len(youngStarIndex))
young_x_pos = []
young_y_pos = []
young_z_pos = []
young_mass = []
young_current_mass = []
young_metals = []
young_age = []
young_x_vel = []
young_y_vel = []
young_z_vel = []
young_smooth = []
young_SFR = []
young_logC = []
young_p = []
young_f_PDR = []
def mass_prob(x): # star cluster mass probability distribution
    return x**(-1.8)
mass_min = 700
mass_max = 1e6
prob_max = 7.57*10**(-6) # slightly larger than 700**(-1.8)
N_MC = len(youngStarIndex)*10000 # number of Monte Carlo samples (yields ~ N_MC/1000 masked samples)
# accept / reject Monte Carlo:
mass = np.random.uniform(mass_min,mass_max,N_MC)  # get uniform temporary mass values
prob = np.random.uniform(0,prob_max,N_MC)  # get uniform random probability values
mask = prob < mass_prob(mass) # accept / reject
sampled_masses = mass[mask] # sample of star cluster masses following the desired distribution (for calculating young_p)
for i in range(len(youngStarIndex)):
    parent_index = youngStarIndex[i]
    parent_mass = cut_mass[parent_index]
    parent_age = cut_age[parent_index]
    young_x_pos.append(cut_x_pos[parent_index])
    young_y_pos.append(cut_y_pos[parent_index])
    young_z_pos.append(cut_z_pos[parent_index])
    young_current_mass.append(cut_mass[parent_index]) # assume no mass loss 
    young_metals.append(cut_metals[parent_index])
    young_age.append(cut_age[parent_index])
    young_x_vel.append(cut_x_vel[parent_index])
    young_y_vel.append(cut_y_vel[parent_index])
    young_z_vel.append(cut_z_vel[parent_index])
    young_smooth.append(cut_smooth[parent_index])
    young_f_PDR.append(1.) 
    age_bins = np.linspace(0.,1.e7,num=1000)
    temp_young_f_PDR = np.trapz(np.exp(-1.*age_bins/tau_clear), x=age_bins) / 1e7
    young_mass.append(parent_mass * temp_young_f_PDR)
    young_SFR.append(1.e-7 * young_mass[i]) # (units: yr**-1) assumes constant SFR over the last 10 Myrs
    young_logC.append(np.random.normal(5., 0.4))
    young_p.append(k*10**((5/2)*(young_logC[i]-(3/5)*np.log10(sampled_masses[i])))) # in kg cm**-1 s**-2
    young_p[i] *= 100 # convert to Pascals
    # set parent mass to 0
    cut_mass[parent_index] = 0.

# create negative mass gas particles to compensate for dust included in MAPPINGS-III SEDs
length = len(young_x_pos)
cut_x_pos_gas = np.append(cut_x_pos_gas, np.asarray(young_x_pos))
cut_y_pos_gas = np.append(cut_y_pos_gas, np.asarray(young_y_pos))
cut_z_pos_gas = np.append(cut_z_pos_gas, np.asarray(young_z_pos))
cut_smooth_gas = np.append(cut_smooth_gas, 3.*np.asarray(young_smooth)) # 3 times larger smoothing length
cut_mass_gas = np.append(cut_mass_gas, -10.*np.asarray(young_mass)) # 10 times larger negative mass
cut_metals_gas = np.append(cut_metals_gas, np.asarray(young_metals))
cut_temp_gas = np.append(cut_temp_gas, np.zeros(length)+8000.) # all have temperature of 8000K (arbitrary)

# save text files
star_header = 'Column 1: position x (pc)\nColumn 2: position y (pc)\nColumn 3: position z (pc)\nColumn 4: smoothing length (pc)\nColumn 5: v_x (km/s)\nColumn 6: v_y (km/s)\nColumn 7: v_z (km/s)\nColumn 8: mass (Msun)\nColumn 9: metallicity ()\nColumn 10: age (yr)'
gas_header = 'Column 1: position x (pc)\nColumn 2: position y (pc)\nColumn 3: position z (pc)\nColumn 4: smoothing length (pc)\nColumn 5: mass (Msun)\nColumn 6: metallicity ()\nColumn 7: temperature (K)' 

np.savetxt(textPath+'stars.txt',np.float32(np.c_[cut_x_pos, cut_y_pos, cut_z_pos, cut_smooth, cut_x_vel, cut_y_vel, cut_z_vel, cut_mass, cut_metals, cut_age]),header=star_header)
np.savetxt(textPath+'gas.txt',np.float32(np.c_[cut_x_pos_gas, cut_y_pos_gas, cut_z_pos_gas, cut_smooth_gas, cut_mass_gas, cut_metals_gas, cut_temp_gas]),header=gas_header)
np.savetxt(textPath+'youngStars.txt',np.float32(np.c_[young_x_pos, young_y_pos, young_z_pos, young_smooth, young_x_vel, young_y_vel, young_z_vel, young_SFR, young_metals, young_logC, young_p, young_f_PDR]))


print('done')




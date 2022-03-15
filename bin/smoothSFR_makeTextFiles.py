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


# create MAPPINGS-III particles
tau_clear = 2.e6 # assume tau_clear years clearing time (e-folding time)
smooth_time = 0. # star-formation lasts smooth_time years
k = 1.3806485*10**(-19) # boltzmann constant in cm**2 kg s**-2 K**-1
#youngStarMask = self.age < 1.1e8 # Mask for young stars (< 110 Myrs)
youngStarMask = cut_age < (smooth_time + 1.e7) # Mask for young stars (smooth_time + 10 Myrs)
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
sampled_masses = mass[mask] # sample of star cluster masses following the desired distribution
for i in range(len(youngStarIndex)):
    parent_index = youngStarIndex[i]
    parent_mass = cut_mass[parent_index]
    parent_age = cut_age[parent_index]
    # parent_mass_not_formed is the fraction of mass that hasn't formed yet
    # parent_mass_fraction is the fraction of mass that has formed that is younger than 10 Myrs
    if (parent_age > 1.e7) and (parent_age <= smooth_time): # 10 Myrs < age < smooth_time
        parent_mass_fraction = 0.1 # 10 percent of the mass is between 0 and 10 Myrs old
        parent_mass_not_formed = (smooth_time - parent_age) / smooth_time
        age_min = 0 
        age_max = 1.e7
    elif parent_age <= 1.e7: # younger than 10 Myrs
        parent_mass_fraction = parent_age / 1.e7
        parent_mass_not_formed = (smooth_time - parent_age) / smooth_time
        age_min = 0
        age_max = parent_age
    else: # between smooth_time Myrs and 110 Myrs
        parent_mass_not_formed = 0.
        age_min = parent_age - smooth_time 
        age_max = 1.e7
        parent_mass_fraction = (age_max - age_min) / 1.e7 
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
    age_bins = np.linspace(age_min,age_max,num=1000)
    temp_young_f_PDR = np.trapz(np.exp(-1.*age_bins/tau_clear), x=age_bins) / 1e6
    temp_young_f_PDR /= (age_max - age_min)/1e6 # average f_PDR over age range in Myrs
    young_mass.append(parent_mass * (1. - parent_mass_not_formed) * parent_mass_fraction * temp_young_f_PDR)
    young_SFR.append(1.e-7 * young_mass[i]) # (units: yr**-1) assumes constant SFR over the last 10 Myrs
    young_logC.append(np.random.normal(5., 0.4))
    young_p.append(k*10**((5/2)*(young_logC[i]-(3/5)*np.log10(sampled_masses[i])))) # in kg cm**-1 s**-2
    young_p[i] *= 100 # convert to Pascals
    # subtract mass from parent particles
    cut_mass[parent_index] -= (cut_mass[parent_index]*parent_mass_not_formed + young_mass[i])

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




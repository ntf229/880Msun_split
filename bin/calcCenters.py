# Calculate and store stellar center of masses for each snapshot

import numpy as np
import os
import matplotlib.pyplot as plt
import h5py

#filePath = '/scratch/ntf229/880Msun_split/uv_background/m12i_res7100_uvb-late/1Myr_res880/output/snapdir_696/'

storePath = '/scratch/ntf229/880Msun_split/resources/centers/'

# Create directories if they don't already exist
os.system('mkdir -p '+storePath)

snaps = ["snapdir_570","snapdir_576","snapdir_581","snapdir_586","snapdir_591","snapdir_596","snapdir_601","snapdir_606","snapdir_611","snapdir_616","snapdir_621","snapdir_626","snapdir_631","snapdir_636","snapdir_641","snapdir_646","snapdir_651","snapdir_656","snapdir_661","snapdir_666","snapdir_671","snapdir_676","snapdir_681","snapdir_686","snapdir_691","snapdir_696"]

boxSize = 6.e4 # in pc

center = np.asarray([4.1825e7, 4.415e7, 4.6275e7]) # center of snapshot 696 in pc (selected by hand)
avgV = np.asarray([0., 0., 0.]) # in km/s

# constants
proton_mass = 8.4089382e-58 # in solar masses
k_Boltzmann = 6.94169e-60 # in km**2 M_sun s**-2 K**-1
gamma = 5./3.
XH = 0.76 
y_helium = (1.0-XH)/(4.0*XH)  # from xiangchang ma physics.py

for s in range(len(snaps)):
    index = (len(snaps) - s - 1) # in reverse order
    currentSnap = snaps[index]
    snap_num = snaps[index].split('_')[1]
    print('starting snapshot', snap_num)
    if index == 0: # parent snapshot
        filePath = '/scratch/ntf229/880Msun_split/uv_background/m12i_res7100_uvb-late/output/'+snaps[index]+'/'
    else:
        filePath = '/scratch/ntf229/880Msun_split/uv_background/m12i_res7100_uvb-late/1Myr_res880/output/'+snaps[index]+'/'
    # move center by assuming constant average stellar velocity since the previous (more recent) snapshot
    avgV *= 3.241e-14 # conver to pc/s
    timeInt = 696. - float(snap_num) # time interval between current snapshot and snapshot 696
    timeInt *= 3.154e13 # convert from Myrs to seconds
    center += -1. * avgV * timeInt # negative because we're going backwards in time
    #np.save(storePath+snap_num+'_center.npy', center)
    
    # Spatial cuts (in pc)
    xAbsMin = center[0] - (boxSize/2)
    xAbsMax = center[0] + (boxSize/2)
    yAbsMin = center[1] - (boxSize/2)
    yAbsMax = center[1] + (boxSize/2)
    zAbsMin = center[2] - (boxSize/2)
    zAbsMax = center[2] + (boxSize/2)
    
    f = h5py.File(filePath+'snapshot_'+snap_num+'.0.hdf5', 'r')
    numFiles = f['Header'].attrs['NumFilesPerSnapshot']
    hubble = f['Header'].attrs['HubbleParam']
    hinv = 1.0 / hubble
    OmegaM = 0.272 
    
    for i in range(numFiles):
        
        # load data from current file
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
    
        # append data
        if i == 0:
            full_x_pos = x_pos 
            full_y_pos = y_pos
            full_z_pos = z_pos
            full_x_vel = x_vel
            full_y_vel = y_vel
            full_z_vel = z_vel
            full_mass = mass
        else:
            full_x_pos = np.append(full_x_pos, x_pos)
            full_y_pos = np.append(full_y_pos, y_pos)
            full_z_pos = np.append(full_z_pos, z_pos)
            full_x_vel = np.append(full_x_vel, x_vel)
            full_y_vel = np.append(full_y_vel, y_vel)
            full_z_vel = np.append(full_z_vel, z_vel)
            full_mass = np.append(full_mass, mass)

    # spatial cuts    
    xMin = full_x_pos > xAbsMin
    xMax = full_x_pos[xMin] < xAbsMax
    yMin = full_y_pos[xMin][xMax] > yAbsMin
    yMax = full_y_pos[xMin][xMax][yMin] < yAbsMax
    zMin = full_z_pos[xMin][xMax][yMin][yMax] > zAbsMin
    zMax = full_z_pos[xMin][xMax][yMin][yMax][zMin] < zAbsMax
    
    cut_x_pos = full_x_pos[xMin][xMax][yMin][yMax][zMin][zMax]
    cut_y_pos = full_y_pos[xMin][xMax][yMin][yMax][zMin][zMax]
    cut_z_pos = full_z_pos[xMin][xMax][yMin][yMax][zMin][zMax]
    cut_x_vel = full_x_vel[xMin][xMax][yMin][yMax][zMin][zMax]
    cut_y_vel = full_y_vel[xMin][xMax][yMin][yMax][zMin][zMax]
    cut_z_vel = full_z_vel[xMin][xMax][yMin][yMax][zMin][zMax]
    cut_mass = full_mass[xMin][xMax][yMin][yMax][zMin][zMax]
    
    avgV = np.asarray([np.average(cut_x_vel), np.average(cut_y_vel), np.average(cut_z_vel)])
    # calculate and store stellar center of mass after spatial cuts
    center = np.asarray([np.average(cut_x_pos, weights=cut_mass), np.average(cut_y_pos, weights=cut_mass), np.average(cut_z_pos, weights=cut_mass)])
    np.save(storePath+snap_num+'.npy', center)

print('done')




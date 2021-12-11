# Calculate and store bandpass fluxes as numpy arrays 

import numpy as np
import pts.simulation as sm
from sedpy import observate
import astropy.io.fits as fits
import os
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb

filePath = '/scratch/ntf229/880Msun_split/SKIRT/'

bands = ["sdss_r0","sdss_z0","sdss_g0"]
#bands = ["wise_w1","wise_w2","wise_w3","wise_w4","galex_NUV"]

os.system('mkdir -p '+filePath+'grids')
wavelengths = np.load(filePath+'wave.npy') * 1e4 # convert from microns to A
hdul = fits.open(filePath+'sph_i00_total.fits')
data = hdul[0].data
print('data shape:',data.shape) # data shape: (601, 2000, 2000)
print('num wavelengths:',data.shape[0])
print('num x pixels:',data.shape[1])
print('num y pixels:',data.shape[2])

filters = observate.FilterSet(bands)

print('max flux:',np.amax(data))

# initialize grid for band magnitudes
grid = np.zeros( (data.shape[1], data.shape[2], len(bands)) ) 

for i in range(data.shape[1]):
    for j in range(data.shape[2]):
        spec = 1e6 * data[:,i,j] # convert from MJy Jy
        print('spec max:',np.amax(spec))
        if np.sum(spec) == 0:
            grid[i,j,:] = 0
        else:
            print('non-zero')
            f_lambda_cgs = (1/33333) * (1/(wavelengths**2)) * (spec) 
            mags = filters.get_sed_maggies(f_lambda_cgs, sourcewave=wavelengths) # maggies
            for b in range(len(bands)):
                grid[i,j,b] = mags[b]
grid[np.isnan(grid)] = 0 # change all nans to 0 
grid[grid == np.inf] = 0 # change all inf to 0
grid[grid == -np.inf] = 0 # change all -inf to 0

for b in range(len(bands)):
    if os.path.isfile(filePath+'grids/'+bands[b]+'.npy'):
        os.system('rm '+filePath+'grids/'+bands[b]+'.npy')
    np.save(filePath+'grids/'+bands[b]+'.npy', grid[:,:,b])
                



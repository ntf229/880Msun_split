import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--index") # index of snapshot (see snaps below)
args = parser.parse_args()
index = int(args.index)

zrg = True # RGB = zrg, else RGB = w4,z,fuv

snaps = ["snapdir_570","snapdir_576","snapdir_581","snapdir_586","snapdir_591","snapdir_596","snapdir_601","snapdir_606","snapdir_611","snapdir_616","snapdir_621","snapdir_626","snapdir_631","snapdir_636","snapdir_641","snapdir_646","snapdir_651","snapdir_656","snapdir_661","snapdir_666","snapdir_671","snapdir_676","snapdir_681","snapdir_686","snapdir_691","snapdir_696"]
snap_num = snaps[index].split('_')[1]
band_names = ['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', '2MASS_J', '2MASS_H', '2MASS_KS', 'W1', 'W2', 'W3', 'W4', 'PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500']

SKIRTpath = '/scratch/ntf229/880Msun_split/SKIRT/'+snaps[index]+'/'
if zrg:
    plotPath = '/scratch/ntf229/880Msun_split/RGB/zrg/'
else:
    plotPath = '/scratch/ntf229/880Msun_split/RGB/w4_z_fuv/'

os.system('mkdir -p '+plotPath)

plt.figure(figsize=(10,8))

file = fits.open(SKIRTpath+'sph_broadband_total.fits') 
data = np.asarray(file[0].data) # (2000, 2000, 20) bands in Jy

os.system('mkdir -p '+plotPath)

if zrg:
    bands = [6,4,3] # RGB = zrg
else:
    bands = [13,6,0]

r_grid = np.asarray(data[bands[0],:,:])
g_grid = np.asarray(data[bands[1],:,:])
b_grid = np.asarray(data[bands[2],:,:])

print('r_grid shape:', r_grid.shape)

if zrg:
    tot_max = np.amax([np.amax(r_grid),np.amax(g_grid),np.amax(b_grid)])
    r_grid = 255 * r_grid / tot_max
    g_grid = 255 * g_grid / tot_max
    b_grid = 255 * b_grid / tot_max
else:
    tot_max = np.amax([np.amax(r_grid),np.amax(g_grid),np.amax(b_grid)])
    r_grid = 255 * r_grid / tot_max
    g_grid = 255 * g_grid / tot_max
    b_grid = 255 * b_grid / tot_max
    # custom scaling
    g_grid = g_grid / 10
    b_grid = b_grid * 20

# Scaled such that average pixel is white
#r_grid = 255 * r_grid / np.amax(r_grid)
#g_grid = 255 * g_grid / np.amax(g_grid)
#b_grid = 255 * b_grid / np.amax(b_grid)
# custom scaling for w4,z,fuv bands
#g_grid = g_grid / 20
#b_grid = b_grid / 4

sizes = np.shape(r_grid)

fig = plt.figure()
fig.set_size_inches(1. * sizes[0] / sizes[1], 1, forward = False)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
if zrg:
    image = make_lupton_rgb(r_grid, g_grid, b_grid, stretch=3, Q=8) # optical stretching
else:
    image = make_lupton_rgb(r_grid, g_grid, b_grid, stretch=0.2, Q=40) # fuv - ir stretching
ax.imshow(image,interpolation='none')
#plt.savefig(plotPath+snap_num+'_'+band_names[bands[0]]+'_'+band_names[bands[1]]+'_'+band_names[bands[2]]+'.png',dpi=300)
plt.savefig(plotPath+snap_num+'_'+band_names[bands[0]]+'_'+band_names[bands[1]]+'_'+band_names[bands[2]]+'.png',dpi=sizes[0])
plt.close()












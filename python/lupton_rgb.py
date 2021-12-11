import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb

plotPath = '/scratch/ntf229/880Msun_split/lupton_images/'
filePath = '/scratch/ntf229/880Msun_split/SKIRT/'
os.system('mkdir -p '+plotPath)

bands = ["sdss_z0","sdss_r0","sdss_g0"] # corresponding to RGB channels
#bands = ["herschel_pacs_100","sdss_g0","galex_FUV"]

r_grid = np.load(filePath+'grids/'+bands[0]+'.npy')
g_grid = np.load(filePath+'grids/'+bands[1]+'.npy')
b_grid = np.load(filePath+'grids/'+bands[2]+'.npy')

print('r grid max:',np.amax(r_grid))
print('g grid max:',np.amax(g_grid))
print('b grid max:',np.amax(b_grid))

tot_max = np.amax([np.amax(r_grid),np.amax(g_grid),np.amax(b_grid)]) 
r_grid = 255 * r_grid / tot_max
g_grid = 255 * g_grid / tot_max
b_grid = 255 * b_grid / tot_max

# Scaled such that average pixel is white
#r_grid = 255 * r_grid / np.amax(r_grid)
#g_grid = 255 * g_grid / np.amax(g_grid)
#b_grid = 255 * b_grid / np.amax(b_grid)

plt.figure(figsize=(2000/300, 2000/300), dpi=300)
image = make_lupton_rgb(r_grid, g_grid, b_grid, stretch=3, Q=8)
plt.axis('off')
plt.imshow(image,interpolation='none')
plt.savefig(plotPath+bands[0]+'_'+bands[1]+'_'+bands[2]+'.png',dpi=300)
plt.close()





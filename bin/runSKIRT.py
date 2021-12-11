import os

SKIRTPath = '/scratch/ntf229/880Msun_split/SKIRT/'
os.system('mkdir -p '+SKIRTPath)
os.chdir(SKIRTPath)
os.system('skirt sph.ski')
os.system('python -m pts.do plot_seds .')
os.system('python -m pts.do make_images .')

print('done')

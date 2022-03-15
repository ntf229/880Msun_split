import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--index") # index of snapshot (see snaps below, 0-25)                                                   $
args = parser.parse_args()
index = int(args.index)

snaps = ["snapdir_570","snapdir_576","snapdir_581","snapdir_586","snapdir_591","snapdir_596","snapdir_601","snapdir_606","snapdir_611","snapdir_616","snapdir_621","snapdir_626","snapdir_631","snapdir_636","snapdir_641","snapdir_646","snapdir_651","snapdir_656","snapdir_661","snapdir_666","snapdir_671","snapdir_676","snapdir_681","snapdir_686","snapdir_691","snapdir_696"]

mainPath = '/home/ntf229/CCA/880Msun_split/'
textPath = '/scratch/ntf229/880Msun_split/TextFiles/'+snaps[index]+'/'
SKIRTPath = '/scratch/ntf229/880Msun_split/SKIRT/'+snaps[index]+'/'

os.system('mkdir -p '+SKIRTPath)
os.chdir(SKIRTPath)
os.system('cp '+textPath+'stars.txt .')
os.system('cp '+textPath+'gas.txt .')
os.system('cp '+textPath+'youngStars.txt .')
os.system('cp /home/ntf229/CCA/880Msun_split/ski_files/sph.ski .')

os.system('python '+mainPath+'python/modify_ski.py --filePath='+SKIRTPath+'/sph.ski --inc=0 --numPhotons=1e9 --pixels=2000')

os.system('skirt sph.ski')
os.system('python -m pts.do plot_seds .')
os.system('python -m pts.do make_images .')
os.system('rm stars.txt')
os.system('rm gas.txt')
os.system('rm youngStars.txt')

print('done')

import yt
import numpy as np
import sys, os

snap_dir = sys.argv[1]
outfile = sys.argv[2]


pos = {}
ngalaxies = {}

#we only want the TNG galaxies that have a stellar mass above the minimum simba stellar mass and that has a nonzero gas mass
#essentially eliminating the differences between the TNG FOF and simba w/caesar
gal_list = np.load('/orange/narayanan/s.lower/TNG/galaxy_list_mass_req.npz')['galaxies']

i = 0
for subhalo in gal_list:
    ds = yt.load(snap_dir+'/snap99_galaxy'+'{:03d}'.format(i)+'.hdf5')
    ad = ds.all_data()
    pos['galaxy'+str(i)] = {}
    xpos = np.sum(ad[('all', 'Coordinates')][:,0].value * ad[('all', 'Masses')].value) / np.sum(ad[('all', 'Masses')].value)
    ypos = np.sum(ad[('all', 'Coordinates')][:,1].value * ad[('all', 'Masses')].value) / np.sum(ad[('all', 'Masses')].value)
    zpos = np.sum(ad[('all', 'Coordinates')][:,2].value * ad[('all', 'Masses')].value) / np.sum(ad[('all', 'Masses')].value)


    pos['galaxy'+str(i)]['snap099'] = np.array([x_pos, y_pos, z_pos])
    i += 1

ngalaxies['snap099'] = len(gal_list)

np.savez(outfile, ngalaxies=ngalaxies, pos=pos)

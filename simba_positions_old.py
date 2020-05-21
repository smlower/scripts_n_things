import yt
import numpy as np
import sys, os
import caesar
import numpy as np
import glob
snap = int(sys.argv[1])



ids = np.load('/orange/narayanan/s.lower/simba/snap'+str(snap)+'_galaxy.npz', allow_pickle=True)
gal = ids['galid']
snap_dir = '/orange/narayanan/s.lower/simba/desika_filtered_snaps/snap'+str(snap)+'/'
outfile = '/orange/narayanan/s.lower/simba/desika_filtered_snaps/simba_snap'+str(snap)+'_pos.npz'
print('loading yt ds')
ds = yt.load('/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_'+"{:03d}".format(snap)+'.hdf5')
ad = ds.all_data()
print('caesar quick loading')
c_files = glob.glob('/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0'+"{:03d}".format(snap)+'_z*.hdf5')
obj = caesar.quick_load(c_files[0])

gas_masses = ad[('PartType0', 'Masses')].value
#star_masses = ad[('PartType4', 'Masses')].value
star_masses = ad[('PartType4', 'Masses')].value
gas_coords = ad[('PartType0', 'Coordinates')].value#[:,0]
star_coords = ad[('PartType4', 'Coordinates')].value#[:,0]


pos = {}
ngalaxies = {}


for i in range(len(gal)):
    print('on galaxy ',i,' of ',len(gal))
    try:
        glist = obj.galaxies[gal[i]].glist
        slist = obj.galaxies[gal[i]].slist
    except:
        print('galaxy',i,'does not exist')
        continue
    pos['galaxy'+str(i)] = {}
    
    total_mass = np.sum(gas_masses[glist]) + np.sum(star_masses[slist])
    
    x_pos = (np.sum(gas_masses[glist] * gas_coords[:,0][glist]) + np.sum(star_masses[slist] * star_coords[:,0][slist])) / total_mass 
    y_pos = (np.sum(gas_masses[glist] * gas_coords[:,1][glist]) + np.sum(star_masses[slist] * star_coords[:,1][slist])) / total_mass
    z_pos = (np.sum(gas_masses[glist] * gas_coords[:,2][glist]) + np.sum(star_masses[slist] * star_coords[:,2][slist])) / total_mass


    pos['galaxy'+str(i)]['snap'+str(snap)] = np.array([x_pos, y_pos, z_pos])
    

ngalaxies['snap'+str(snap)] = len(obj.galaxies)

np.savez(outfile, ngalaxies=ngalaxies, pos=pos)

from pyloser.los_extinction import set_solar, compute_AV, init_kerntab
import yt
import sys
import numpy as np

galaxy = sys.argv[1]


galaxy_file = '/orange/narayanan/s.lower/simba/filtered_snapshots/snap305/galaxy_'+str(galaxy)+'.hdf5'

#load galaxy snapshot
print('loading snapshot data')
ds = yt.load(galaxy_file)
ad = ds.all_data()
h = ds.hubble_constant
#load particle info
print('loading particle info')
gas = 'PartType0'
stars = 'PartType4'
smass = ad[(stars, 'Masses')]
gmass = ad[(gas,'Masses')].in_units('Msun').value 
gpos = ad[(gas, 'Coordinates')].in_units('kpc').value
ghsm = ad[(gas, 'SmoothingLength')].in_units('kpc').value
spos = ad[(stars, 'Coordinates')].in_units('kpc').value
#dust_masses = ad[(gas, 'Dust_Masses')].in_units('Msun')
dustmass=ds.arr(ad[('PartType0','Dust_Masses')],'code_mass')

gmet = dustmass.in_units('Msun').value / gmass



#setup minimum shell param file for load_sim to read
params = {}
params['map_type'] = 'object'
params['output_type'] = 'phot'
params['gas_metal_index'] = -1


#setup aux variables for los extinction
nproc = 1
viewing_direction = 0
kerntab = init_kerntab('cubic')
Solar = set_solar()
lowZ_scaling = 'Simba_DTM'

    

print('computing AV')
AV_per_star = compute_AV(viewing_direction,gmass,gpos,ghsm,gmet,spos,Solar,25000., ds.scale_factor, 'cubic',kerntab,lowZ_scaling)
    
print('saving')
np.savez('/orange/narayanan/s.lower/simba/pyloser_los_extinction/pyloser_computeAV_galaxy'+str(galaxy)+'.npz', AV_per_star=AV_per_star, stars=len(smass))





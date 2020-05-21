import sys
import numpy as np
import yt


#ids = np.load('/orange/narayanan/s.lower/simba/snap'+str(args.snapnum)+'_galaxy.npz', allow_pickle=True)
#gal = ids['galid'][args.galaxy]

galaxy = sys.argv[1]
print('loading snapshot')
filt_snaps = '/orange/narayanan/s.lower/simba/filtered_snapshots/snap305/'
ds = yt.load(filt_snaps+'/galaxy_'+str(galaxy)+".hdf5")
data = ds.all_data()
print('getting dust mass')
dust_mass = data.ds.arr(data[("PartType0", "Dust_Masses")].value, 'code_mass')
true_dust_mass = np.sum(dust_mass.in_units('Msun').value)


np.savez(filt_snaps+'/galaxy'+str(galaxy)+'_dustmass.npz', galaxy=galaxy, dust_mass=true_dust_mass)

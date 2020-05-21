import yt
import numpy as np
import readhaloHDF5 as readhalo


h = readhalo.HaloReader('/orange/narayanan/s.lower/', '099', 99)


simba_mstar_lim = 4.4e7 #Msun

subhalo_idx = []

for subhalo in range(len(ad['Subhalo', 'particle_identifier'])):
    star_mass = ad['Subhalo', 'SubhaloMassType_4'][subhalo].in_units('Msun')
    gas_mass = ad['Subhalo', 'SubhaloMassType_0'][subhalo].in_units('Msun')
    if (star_mass.value > simba_mstar_lim) and (gas_mass.value > 0):
        subhalo_idx.append(ad['Subhalo', 'particle_identifier'][subhalo])
    

np.savez('/orange/narayanan/s.lower/output/tng_mass_cutoff.npz', galaxy=subhalo_idx)

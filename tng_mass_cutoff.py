import yt
import numpy as np
import sys, os
import readhaloHDF5 as readhalo


h = readhalo.HaloReader('/orange/narayanan/s.lower/TNG/', '099', 99)

#fof_dir = sys.argv[1]

#ds = yt.load(fof_dir+'fof_subhalo_tab_099.0.hdf5')

#ad = ds.all_data()

simba_mstar_lim = 3.0e8 #Msun

subhalo_idx = []

#for subhalo in range(len(ad['Subhalo', 'particle_identifier'])):
#    star_mass = ad['Subhalo', 'SubhaloMassType_4'][subhalo].in_units('Msun')
#    gas_mass = ad['Subhalo', 'SubhaloMassType_0'][subhalo].in_units('Msun')
#    if (star_mass.value > simba_mstar_lim) and (gas_mass.value > 0):
#        subhalo_idx.append(ad['Subhalo', 'particle_identifier'][subhalo])
    

#print(len(subhalo_idx))



from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
#import pdb,ipdb

import pickle

import sys
from decimal import Decimal

#sys.path.insert(0,'/home/desika.narayanan/zooms/gizmo_stuff/analysis/caesar/')


file = '/ufrc/narayanan/s.lower/simSEDs/m25_n512_simba_prospector/physical_properties/0_499.Snapshot305.pkl'

with open(file,'rb') as pickle_file: content = pickle.load(pickle_file)

color_array = ['navy','cornflowerblue','mediumseagreen']


fig = plt.figure()
ax = fig.add_subplot(111)

#first plot for a few random halos

#for counter,halo_idx in enumerate([0, 10, 95]):

    #z = content.halos[halo_idx].z
    #central_sfh = content.halos[halo_idx].central_sfh
    #h_sfh = content.halos[halo_idx].h_sfh

    #get the mstar and then make it nice format
    #z0_mstar = content.halos[halo_idx].central_mstar[1] #we use the 1st index since the 0th index is actually not saved
    #z0_mstar_label = '%.1E' % Decimal(float(z0_mstar.value))

    #ax.plot(z,central_sfh,label=r'z=0 M$_*$ = '+z0_mstar_label,color=color_array[counter])
    #ax.plot(z,h_sfh,label='Halo SFR')

#ax.set_yscale('log')
#ax.set_ylabel(r'SFR (M$_\odot$ yr$^{-1}$')
#ax.set_xlabel('z')
#ax.set_xlim(ax.get_xlim()[::-1])
#ax.set_ylim([1,1.e4])
#ax.set_xlim([7,0])
#plt.legend(loc=2)
#fig.savefig('ridic_sfr.png',dpi=300)


#now just for halo0 plot the halo SFR and the central sfr

for counter,halo_idx in enumerate(np.arange(10)):

    z = content.halos[halo_idx].z
    central_sfh = content.halos[halo_idx].central_sfh
    h_sfh = content.halos[halo_idx].h_sfh

    print(central_sfh)


    #get the mstar and then make it nice format
    #z0_mstar = content.halos[halo_idx].central_mstar[0] #we use the 1st index since the 0th index is actually not saved
    #z0_mstar_label = '%.1E' % Decimal(float(z0_mstar.value))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(z,central_sfh,label=r'z=0 halo '+str(halo_idx),color=color_array[0])
    #ax.plot(z,h_sfh,'--',label=r'z=0 halo'+halo_idx,color=color_array[0])
    #ax.fill_between(z,central_sfh,h_sfh,color='orange')

    
    ax.set_yscale('log')
    ax.set_ylabel(r'SFR (M$_\odot$ yr$^{-1}$')
    ax.set_xlabel('z')
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylim([1,1.e4])
    ax.set_xlim([7,0])
    plt.legend(loc=2)
    fig.savefig('sfh_halos_centrals'+str(halo_idx)+'.png',dpi=300)

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import caesar
#import pdb,ipdb

from astropy import units as u
from astropy import constants as const
from glob2 import glob

from get_pandas_physical_quantities import *

import yt
import pickle

#-----------------------------------
#MODIFIABLE HEADER

snapshot_directory = '/ufrc/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/'
g_outfile = 'm25n512.galaxies.pandas.hdf5'
h_outfile = 'm25n512.halos.pandas.hdf5'


NMAXHALOS = 1.e4 #we only save the top 10,000 halos from the lowest redshift snapshot
CHUNKSIZE = 500 #number of galaxies per file
#-----------------------------------
CHUNKSIZE = int(CHUNKSIZE)

def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)



def gadget_snap_num(snapnum):
    if snapnum < 10:
        gadget_snap_num = '00'+str(snapnum)
    elif snapnum >= 10 and snapnum <100:
        gadget_snap_num = '0'+str(snapnum)
    else:
        gadget_snap_num = str(snapnum)
    return gadget_snap_num


MEMBERS = np.sort(glob('%s/caesar*.hdf5' % (snapshot_directory+'/Groups_oldprogen/')))
try:  LASTSNAP = int(MEMBERS[-1][-16:-12])
except ValueError:LASTSNAP = int(MEMBERS[-1][-8:-5])



g_data = {}
h_data = {}
snaplist = []

for fileidx, file in enumerate(MEMBERS[::-1][0:270]):

    print('loading %s' % file)

    obj = caesar.load(file)
    snapshotname = obj.simulation.basename
    g_mstar = [gal.masses['stellar'] for gal in obj.galaxies]
    g_sfr = [gal.sfr for gal in obj.galaxies]
    redshift = obj.simulation.redshift
    g_mgas = [gal.masses['gas'] for gal in obj.galaxies]
    g_mgas_H2 = [gal.masses['H2'] for gal in obj.galaxies]
    g_mgas_HI = [gal.masses['HI'] for gal in obj.galaxies]
    g_mdust = [gal.masses['dust'] for gal in obj.galaxies]

    try: g_progen_index = [gal.progen_index for gal in obj.galaxies]
    except AttributeError: g_progen_index = -1
    parent_halo_index = [gal.parent_halo_index for gal in obj.galaxies]

    h_mstar = [halo.masses['stellar'] for halo in obj.halos]
    h_mdm = [halo.masses['dm'] for halo in obj.halos]
    h_sfr = [halo.sfr for halo in obj.halos]
    h_mgas = [halo.masses['gas'] for halo in obj.halos]
    h_mgas_H2 = [halo.masses['H2'] for halo in obj.halos]
    h_mgas_HI = [halo.masses['HI'] for halo in obj.halos]
    h_mdust = [halo.masses['dust'] for halo in obj.halos]

    try:h_progen_index = [halo.progen_index for halo in obj.halos]
    except AttributeError: h_progen_index = -1
    galaxy_index_list = [halo.galaxy_index_list for halo in obj.halos]
    
    
    #ew have to do list appening and try/excepts here since not every halo has a central galaxy attribute
    h_central_mstar = []
    h_central_sfr = []
    h_central_mgas = []
    h_central_mgas_H2 = []
    h_central_mgas_HI = []
    h_central_mdust = []
    
    for i in range(obj.nhalos):
        try: 
            h_central_mstar.append(obj.halos[i].central_galaxy.masses['stellar'])
            h_central_mgas.append(obj.halos[i].central_galaxy.masses['gas'])
            h_central_mgas_H2.append(obj.halos[i].central_galaxy.masses['H2'])
            h_central_mgas_HI.append(obj.halos[i].central_galaxy.masses['HI'])
            h_central_sfr.append(obj.halos[i].central_galaxy.sfr)
            h_central_mdust.append(obj.halos[i].central_galaxy.masses['dust'])
        except AttributeError:
            h_central_mstar.append(-1)
            h_central_mgas.append(-1)
            h_central_mgas_H2.append(-1)
            h_central_mgas_HI.append(-1)
            h_central_sfr.append(-1)
            h_central_mdust.append(-1)

    g = {'sfr':g_sfr,'mstar':g_mstar,'redshift':redshift,'mgas':g_mgas,'mgas_H2':g_mgas_H2,'mgas_HI':g_mgas_HI,'mdust':g_mdust,'progen_index':g_progen_index,'parent_halo_index':parent_halo_index}

    h = {'h_sfr':h_sfr,'h_mstar':h_mstar,'h_mdm':h_mdm,'h_mgas':h_mgas,'h_mgas_H2':h_mgas_H2,'h_mgas_HI':h_mgas_HI,'h_mdust':h_mdust,'progen_index':h_progen_index,'galaxy_index_list':galaxy_index_list,'redshift':redshift,'h_central_mstar':h_central_mstar,'h_central_sfr':h_central_sfr,'h_central_mgas':h_central_mgas,'h_central_mgas_H2':h_central_mgas_H2,'h_central_mgas_HI':h_central_mgas_HI,'h_central_mdust':h_central_mdust}

    
    galnames = ['gal'+str(i) for i in range(obj.ngalaxies)]
    g1 = pd.DataFrame(g,index=galnames)

    halonames = ['halo'+str(i) for i in range(obj.nhalos)]
    h1 = pd.DataFrame(h,index=halonames)

    snapnum = int(file[file.find('caesar_')+7:file.find('caesar_')+11])
    snaplist.append(snapnum)

    g_data['Snapshot'+gadget_snap_num(snapnum)] = g1
    h_data['Snapshot'+gadget_snap_num(snapnum)] = h1
    
    

    
g_p1 = pd.Panel(g_data)
h_p1 = pd.Panel(h_data)


#get the physical properties evolution

#first take the minimum of NMAXHALOS and the number of actual halos, just in case there aren't actually NMAXHALOS in the lowest redshift snapshot
NMAXHALOS = np.min([NMAXHALOS,len(h_p1.axes[1])])
halolist = np.arange(NMAXHALOS)
#https://www.geeksforgeeks.org/break-list-chunks-size-n-python/
chunked_list =  [halolist[i * CHUNKSIZE:(i + 1) * CHUNKSIZE] for i in range((len(halolist) + CHUNKSIZE - 1) // CHUNKSIZE )]

snapshots = list(h_p1.axes[0])

for snapshot in snapshots:
    snapnum = int(snapshot[snapshot.find('t')+1:snapshot.find('t')+4])


    for sublist in chunked_list: #iterating through the sublists now
    
        
        halo_list = []
        
        for halonum in sublist:

            z,central_sfh,central_mstar,central_mgas,central_mgas_H2,central_mdust,h_sfh,h_mstar,h_mgas,h_mdm,h_mdust,progen_idx_list = get_quantities(h_p1,snapnum,int(halonum))
            
            halo = haloclass(0,0,0,0,0,0,0,0,0,0,0,0,0)
            halo.z = z
            halo.halonum = halonum
            halo.central_sfh = central_sfh
            halo.central_mstar = central_mstar
            halo.central_mgas = central_mgas
            halo.central_mgas_H2 = central_mgas_H2
            halo.central_mdust = central_mdust
            halo.h_sfh = h_sfh
            halo.h_mstar = h_mstar
            halo.h_mgas = h_mgas
            halo.h_mdm = h_mdm
            halo.h_mdust = h_mdust
            halo.h_progen_idx_list = progen_idx_list
            
            halo_list.append(halo)
            
            masterobj = masterclass(halo_list)
            
            #setting the outfile
        firstgal = int(np.min(sublist))
        lastgal = int(np.max(sublist))
        outfile = '/ufrc/narayanan/s.lower/simSEDs/m25_n512_simba_prospector/physical_properties/'+str(firstgal)+'_'+str(lastgal)+'.'+snapshot+'.pkl'

        save_object(masterobj,outfile)
 
    #read in as such: with open('dum.pkl','rb') as pickle_file: content = pickle.load(pickle_file)



#g_p1.to_hdf(g_outfile,key='df')
#h_p1.to_hdf(h_outfile,key='df')

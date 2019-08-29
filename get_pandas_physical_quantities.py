import pandas as pd
import numpy as np

#import pdb,ipdb


class masterclass:
    def __init__(self,halos):

        self.halos = halos

        
class haloclass:
    def __init__(self,halonum,z,central_sfh,central_mstar,central_mgas,central_mgas_H2,central_mdust,h_sfh,h_mstar,h_mgas,h_mdm,h_mdust,h_progen_idx_list):
        
        self.halonum = 0
        self.z = []
        self.central_sfh = []
        self.central_mstar = []
        self.central_mgas = []
        self.central_mgas_H2 = []
        self.central_dust = []
        self.h_sfh = []
        self.h_mstar = []
        self.h_mgas = []
        self.h_mdm = []
        self.h_mdust = []
        self.h_progen_idx_list = []

def gadget_snap_num(snapnum):
    if snapnum < 10:
        gadget_snap_num = '00'+str(snapnum)
    elif snapnum >= 10 and snapnum <100:
        gadget_snap_num = '0'+str(snapnum)
    else:
        gadget_snap_num = str(snapnum)
    return gadget_snap_num


def link_central_galaxy_snapshots(df,snapnum,halonum):
    #the methodology here is to link the central galaxies from halos
    #and return a list of the indices for prior snapshots
    snapshots = list(df.axes[0])
    

    try: snap_idx = snapshots.index('Snapshot'+gadget_snap_num(snapnum))
    except: print("in [get_pandas_physical_quantities/link_central_galaxy_snapshots the snapshot does not exist in the pandas datafram.   These are the snapshots that exist:\n\n",df.axes[0])

    snapshots_to_link = snapshots[snap_idx:-1]

    progen_idx_list = []
    for snapshot in snapshots_to_link:
        try: progen_idx_list.append(int(df[snapshot]['progen_index']['halo'+str(halonum)]))
        except ValueError:
            progen_idx_list.append(-1)

    return progen_idx_list




def get_quantities(df,snapnum,halonum):
    
    central_sfh = []
    central_mstar = []
    central_mgas = []
    central_mgas_H2 = []
    central_mdust = []
    h_mstar = []
    h_mgas = []
    h_sfh = []
    h_mdm = []
    h_mdust = []
    z = []
    

    snapshots = list(df.axes[0])
    try:snap_idx = snapshots.index('Snapshot'+gadget_snap_num(snapnum))
    except: 
        print("in [get_pandas_physical_quantities/link_central_galaxy_snapshots the snapshot does not exist in the pandas datafram.   These are the snapshots that exist:\n\n",df.axes[0])
        return

    snapshots_to_examine =snapshots[snap_idx:-1]
    progen_idx_list = link_central_galaxy_snapshots(df,snapnum,halonum)

    for counter,snapshot in enumerate(snapshots_to_examine):
        #catch in case we're in a situation where there's no progenitors anymore
        if (halonum == -1) or (progen_idx_list[counter-1] == -1):
            central_sfh.append(-1)
            central_mstar.append(-1)
            central_mgas.append(-1)
            central_mgas_H2.append(-1)
            central_mdust.append(-1)
            h_sfh.append(-1)
            h_mstar.append(-1)
            h_mgas.append(-1)
            h_mdm.append(-1)
            h_mdust.append(-1)
            z.append(-1)
        
        else:
            if counter == 0:
                
                central_sfh.append(df[snapshot]['h_central_sfr']['halo'+str(halonum)])
                central_mstar.append(df[snapshot]['h_central_mstar']['halo'+str(halonum)])
                central_mgas.append(df[snapshot]['h_central_mgas']['halo'+str(halonum)])
                central_mgas_H2.append(df[snapshot]['h_central_mgas_H2']['halo'+str(halonum)])
                central_mdust.append(df[snapshot]['h_central_mdust']['halo'+str(halonum)])
                h_sfh.append(df[snapshot]['h_sfr']['halo'+str(halonum)])
                h_mstar.append(df[snapshot]['h_mstar']['halo'+str(halonum)])
                h_mgas.append(df[snapshot]['h_mgas']['halo'+str(halonum)])
                h_mdm.append(df[snapshot]['h_mdm']['halo'+str(halonum)])
                h_mdust.append(df[snapshot]['h_mdust']['halo'+str(halonum)])
                z.append(df[snapshot]['redshift']['halo'+str(halonum)])
                
            else:
                central_sfh.append(df[snapshot]['h_central_sfr']['halo'+str(progen_idx_list[counter-1])])
                central_mstar.append(df[snapshot]['h_central_mstar']['halo'+str(progen_idx_list[counter-1])])
                central_mgas.append(df[snapshot]['h_central_mgas']['halo'+str(halonum)])
                central_mgas_H2.append(df[snapshot]['h_central_mgas_H2']['halo'+str(halonum)])
                central_mdust.append(df[snapshot]['h_central_mdust']['halo'+str(halonum)])
                h_sfh.append(df[snapshot]['h_sfr']['halo'+str(progen_idx_list[counter-1])])
                h_mstar.append(df[snapshot]['h_mstar']['halo'+str(progen_idx_list[counter-1])])
                h_mgas.append(df[snapshot]['h_mgas']['halo'+str(progen_idx_list[counter-1])])
                h_mdm.append(df[snapshot]['h_mdm']['halo'+str(progen_idx_list[counter-1])])
                h_mdust.append(df[snapshot]['h_mdust']['halo'+str(progen_idx_list[counter-1])])
                z.append(df[snapshot]['redshift']['halo'+str(progen_idx_list[counter-1])])


    return z,central_sfh,central_mstar,central_mgas,central_mgas_H2,central_mdust,h_sfh,h_mstar,h_mgas,h_mdm,h_mdust,progen_idx_list


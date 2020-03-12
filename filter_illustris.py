import h5py
import argparse
import readhaloHDF5 as readhalo
import numpy as np
from glob import glob




parser = argparse.ArgumentParser()
parser.add_argument('base_path')
parser.add_argument('snap_num')
parser.add_argument('galaxy')
parser.add_argument('output')
args = parser.parse_args()



gal_list = np.load('/orange/narayanan/s.lower/TNG/galaxy_list_mass_req.npz')['galaxies']
gal = gal_list[int(args.galaxy)]


h = readhalo.HaloReader(str(args.base_path), '0'+str(args.snap_num), int(args.snap_num))


gas_offset = h.group_offset[gal][0]
star_offset = h.group_offset[gal][4]
gas_len = h.cat.SubhaloLenType[gal][0]
star_len = h.cat.SubhaloLenType[gal][4]



files_dir = args.base_path+'/output/snapdir_'+'0'+str(args.snap_num)+'/'
files = sorted(glob(files_dir+'*.hdf5'), key = lambda name: int(name.split('.')[-2]))
this_file_start0 = 0
this_file_start4 = 0

print('starting with gas particles')
for snap_file in files:
    with h5py.File(snap_file, 'r') as input_file:

        print(snap_file)

        this_file_end0 = this_file_start0 + len(input_file['PartType0']['Masses'])

        #gas_checks
        if (gas_offset+gas_len) > this_file_end0:
            this_file_start0 += len(input_file['PartType0']['Masses'])
            continue

        elif this_file_start0 > (gas_offset+gas_len):
            
            break
        
        else:
            output_file = h5py.File(args.output+'/snap099_galaxy'+"{:03d}".format(int(args.galaxy))+'.hdf5', 'w')
            output_file.copy(input_file['Header'], 'Header')
            output_file.copy(input_file['Config'], 'Config')
            output_file.create_group('PartType0')
            ds_in = yt.load(snap_file)

            print('writing from file',snap_file)
            
            for k in input_file['PartType0']:
                if k in output_file['PartType0']:
                    print('we are writing for the second time, i.e. galaxy data exists across two snapshots')
                    in_data = input_file['PartType0'][k][0 : leftover]
                    output_file['PartType0'][k] = np.hstack(output_file['PartType0'][k], in_data)
                    print(len(in_data))
                else: 
                    print('writing for the first time')
                    in_data = input_file['PartType0'][k][gas_offset-this_file_start0 : min(this_file_end0,(gas_offset+gas_len)-this_file_start0)]
                    output_file['PartType0'][k] = in_data
                    leftover = (gas_offset+gas_len) - this_file_end0
            this_file_start0 += len(input_file['PartType0']['Masses'])

            output_file.close()

print('now doing star particles')
for snap_file in files:
    with h5py.File(snap_file, 'r') as input_file:
        print(snap_file)
        this_file_end4 = this_file_start4 + len(input_file['PartType4']['Masses'])

        #star_checks                                                                           
        if (star_offset+star_len) > this_file_end4:
            this_file_start4 += len(input_file['PartType4']['Masses'])
            continue

        elif this_file_start4 > (star_offset+star_len):
            
            break

        else: 
            output_file = h5py.File(args.output+'/snap099_galaxy'+"{:03d}".format(int(args.galaxy))+'.hdf5')
            output_file.create_group('PartType4')
            print('writing from file',snap_file)
            for k in input_file['PartType4']:
                if k in output_file['PartType4']:
                    print('we are writing for the second time, i.e. galaxy data exists across two files')
                    in_data = input_file['PartType4'][k][0 : leftover]
                    output_file['PartType4'][k] = np.hstack(output_file['PartType4'][k], in_data)
                    
                else:
                    print('writing for the first time')
                    in_data = input_file['PartType4'][k][star_offset-this_file_start4 : min(this_file_end4,(star_offset+star_len)-this_file_start4)]
                    output_file['PartType4'][k] = in_data
                    leftover = (star_offset+gas_len) - this_file_end4
            this_file_start4 += len(input_file['PartType4']['Masses'])


re_out = h5py.File(args.output+'/snap099_galaxy'+"{:03d}".format(int(args.galaxy))+'.hdf5')

re_out['Header'].attrs.modify('NumFilesPerSnapshot', 1)
re_out['Header'].attrs.modify('NumPart_ThisFile', np.array([gas_len, 0, 0, 0, star_len, 0]))
re_out['Header'].attrs.modify('NumPart_Total', np.array([gas_len, 0, 0, 0, star_len, 0]))

print('galaxy '+args.galaxy+' is done.')

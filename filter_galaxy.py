import h5py
import caesar
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ds')
parser.add_argument('caesar')
parser.add_argument('galaxy', type=int)
parser.add_argument('output')
args = parser.parse_args()

obj = caesar.load(args.caesar, LoadHalo=0)

glist = obj.galaxies[args.galaxy].glist
slist = obj.galaxies[args.galaxy].slist
#bhlist = obj.galaxies[args.galaxy].bhlist

with h5py.File(args.ds, 'r') as input_file, h5py.File(args.output, 'w') as output_file:
    output_file.copy(input_file['Header'], 'Header')

    output_file.create_group('PartType0')
    for k in input_file['PartType0']:
        output_file['PartType0'][k] = input_file['PartType0'][k][:][glist]

    output_file.create_group('PartType4')
    for k in input_file['PartType4']:
        output_file['PartType4'][k] = input_file['PartType4'][k][:][slist]

    #output_file.create_group('PartType5')
    #for k in input_file['PartType5']:
    #    output_file['PartType5'][k] = input_file['PartType5'][k][:][bhlist]

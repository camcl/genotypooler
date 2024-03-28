import os
import argparse

inmap = '/home/camille/MagicWheat/src/genotypooler/examples/1_interpolated_wheat_map'
outmap = '/home/camille/MagicWheat/src/genotypooler/examples/1_interpolated_wheat_map_plink.map'
idchrom = 1

linout = []
with open(inmap, 'r') as infile:
    id_pos_cm = infile.readlines()

for ipc in id_pos_cm:
    fields = ipc.split()
    # PLINK format:
    # <Chromosome code> <Variant identifier> <position in cM> <base-pair coordinate>
    linout.append(f'{idchrom} {fields[0]} {fields[2]} {fields[1]}')

with open(outmap, 'w') as outfile:
    for lin in linout:
        outfile.write(lin + '\n')

assert os.path.exists(outmap)

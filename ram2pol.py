# ram2pol.py

# Ilseung Han
#   29.05.2025
#   04.07.2025
#   15.08.2025
#   20.08.2025
#   23.08.2025
#   31.08.2025
#   30.10.2025
#   31.10.2025
#   14.12.2025
#   31.01.2026
#   11.03.2026

import sys
import math
import pymses
import struct
import numpy as np
import astropy.units as u
import astropy.constants as co

from pymses.utils import constants as C
from pymses.filters import CellsToPoints

# constants
grid_id = 20       # grid ID (20 = octree)  
# gamma   = 1.4      # Ramses default
# mH      = 1.66e-24 # Hydrogen mass in g
# kB      = 1.38e-16 # Boltzmann constant in cgs

# X = 0.76 # H mass fraction
# Y = 0.24 # He mass fraction

CLR_LINE =   "                                                      \r"

max_level    = 0
cell_counter = 0
nr_of_cells  = 0

class cell_oct:
    """written by Stefan Reissl"""
    def __init__(self, _x_min, _y_min, _z_min, _length, _level):
        self.x_min = _x_min
        self.y_min = _y_min
        self.z_min = _z_min
        
        self.length = _length
        self.level  = _level
    
        self.isleaf   = 0
        self.data     = []
        self.branches = []  
        
class OcTree:
    """written by Stefan Reissl"""
    def __init__(self, _x_min, _y_min, _z_min, _length):
        self.root = cell_oct(_x_min, _y_min, _z_min, _length, 0)

    def initCellBoundaries(self, cell, _level):
        x_min = cell.x_min
        y_min = cell.y_min
        z_min = cell.z_min
        l     = 0.5 * cell.length

        level = _level

        cell.isleaf      = 0
        cell.data        = []
        cell.branches    = [None, None, None, None, None, None, None, None]
        cell.branches[0] = cell_oct(x_min, y_min, z_min, l, level)
        cell.branches[1] = cell_oct(x_min + l, y_min, z_min, l, level)
        cell.branches[2] = cell_oct(x_min, y_min + l, z_min, l, level)
        cell.branches[3] = cell_oct(x_min + l, y_min + l, z_min, l, level)

        cell.branches[4] = cell_oct(x_min, y_min, z_min + l, l, level)
        cell.branches[5] = cell_oct(x_min + l, y_min, z_min + l, l, level)
        cell.branches[6] = cell_oct(x_min, y_min + l, z_min + l, l, level)
        cell.branches[7] = cell_oct(x_min + l, y_min + l, z_min + l, l, level)     
        
    def insertInTree(self, cell_pos, cell, _level):    
        x_pos = cell.x_min
        y_pos = cell.y_min
        z_pos = cell.z_min

        if cell_pos.level == cell.level:
            cell_pos.data   = cell.data
            cell_pos.isleaf = 1
                        
        else:    
            if len(cell_pos.branches) == 0:
                self.initCellBoundaries(cell_pos, _level + 1)

            x_mid = cell_pos.x_min + 0.5 * cell_pos.length
            y_mid = cell_pos.y_min + 0.5 * cell_pos.length
            z_mid = cell_pos.z_min + 0.5 * cell_pos.length
            
            new_cell_pos = cell_pos

            if(z_pos < z_mid):          # z 0 1 2 3

                if(y_pos < y_mid):      # y 0 1

                    if(x_pos < x_mid):  # x 0
                        new_cell_pos = cell_pos.branches[0]
                    else:               # x 1
                        new_cell_pos = cell_pos.branches[1]

                else:                   # y 2 3

                    if(x_pos < x_mid):  # x 2
                        new_cell_pos = cell_pos.branches[2]
                    else:               # x 3
                        new_cell_pos = cell_pos.branches[3]

            else:                       # z 4 5 6 7

                if(y_pos < y_mid):      # y 4 5

                    if(x_pos < x_mid):  # x 4
                        new_cell_pos = cell_pos.branches[4]
                    else:               # x 5
                        new_cell_pos = cell_pos.branches[5]

                else:                   # y 6 7

                    if(x_pos < x_mid):  # x 6
                        new_cell_pos = cell_pos.branches[6]
                    else:               # x 7
                        new_cell_pos = cell_pos.branches[7]

            self.insertInTree(new_cell_pos, cell, _level + 1)

    def writeOcTree(self, file, cell):
        global cell_counter
        global nr_of_cells

        file.write(struct.pack("H", cell.isleaf))
        file.write(struct.pack("H", cell.level))   

        if cell.isleaf == 1:    
            data_len = len(cell.data)
            
            if cell_counter % 10000 == 0:
                sys.stdout.write('-> Writing octree grid file : ' + str(100.0 * cell_counter / nr_of_cells) + ' %\r')
                sys.stdout.flush()
                
            cell_counter += 1 
         
            for i in range(0, data_len):
                file.write(struct.pack('f', cell.data[i]))
        else:
            for i in range(8):
                self.writeOcTree(file, cell.branches[i])
                
    def checkOcTree(self, cell):
        global cell_counter
        global nr_of_cells

        if cell.isleaf == 1:    
            length = len(cell.data)
            
            if length == 0:
                return False
            
            if cell_counter % 10000 == 0:
                sys.stdout.write('-> Checking octree integrity : ' + str(100.0 * cell_counter / nr_of_cells) + ' %\r')
                sys.stdout.flush()
                
            cell_counter += 1    
            
        else:
            length = len(cell.branches)
            
            if length == 0:
                return False
            
            for i in range(8):
                self.checkOcTree(cell.branches[i])                
                
        return True
        
def loadRamsesData(filename, has_dust_in_sim = True):
    """written by Valeska Valdivia"""
    # filename will be the format my/folder/here/output_?????
    # where ????? is from 00001 to 99999
    # print(filename)
    
    # split the filename into folder and number (pymses needs this)
    outstr = 'output_'
    outloc = filename.rfind(outstr)
    folder = filename[:outloc]
    numloc = outloc + len(outstr)
    
    # note: 5 characters in output number
    num = int(filename[numloc:numloc+5])

    # create the pymses RamsesOutput object 
    # print('\n\nfolder: ', folder , num, '\n\n')
    snap = pymses.RamsesOutput(folder,num)

    # spatial information
    output = {}

    if has_dust_in_sim:
        print("\n--- Loading data for a simulation WITH dust fields ---")
        # load dust-to-fluid mass ratios for all available dust species
        # print(snap.info)
        n_dust = snap.info['ndust']
        dust_ratios = [f'dust_ratio_{i + 1}' for i in range(n_dust)]
    
        # create a flat structure with the snapshot's cell data in it
        amr = snap.amr_source(['rho', 'P', 'vel', 'Br', 'Bl'] + dust_ratios)
        cell_source = CellsToPoints(amr)
        cells = cell_source.flatten()
        # print('The B fields', cells['Br'], cells['Bl'], 'end')
  
        # Calculate densities based on dust ratios. 
        # mass density of fluid
        output['denstot'] = cells['rho'] * snap.info['unit_density'].express(C.g_cc)    # cgs
        output['denstot'] = output['denstot'] * 1e+03                                   # SI
        # print('dens min, max', output['denstot'].min(), output['denstot'].max())

        # mass density of dust and gas
        output['epsilon'] = sum(cells[dust_ratio] for dust_ratio in dust_ratios)
        output['densdto'] = output['denstot'] * output['epsilon']
        for i in range(n_dust):
            output[f'densd{(i + 1):02d}'] = output['denstot'] * cells[f'dust_ratio_{i + 1}']
        output['densgas'] = output['denstot'] - output['densdto']

    else:
        print("\n--- Loading data for a simulation WITHOUT dust fields ---")
        # Assume a single, virtual dust species.
        n_dust = 1
        amr = snap.amr_source(['rho', 'P', 'vel', 'Br', 'Bl'])
        cells = CellsToPoints(amr).flatten()

        # Calculate gas density directly and dust density as 1% of gas density.
        output['densgas'] = cells['rho'] * snap.info['unit_density'].express(C.g_cc)    # cgs
        output['densgas'] = output['densgas'] * 1e+03                                   # SI
        output['densd01'] = output['densgas'] * 1e-02

    # print("\n--- Loading data for a simulation WITH dust fields ---")
    # n_dust = 40
    # amr = snap.amr_source(['rho', 'P', 'vel', 'Br', 'Bl'] +
    #                       [f'rhod{(i + 1):02d}' for i in range(n_dust)])
    # cell_source = CellsToPoints(amr)
    # cells = cell_source.flatten()

    # output['densgas'] = cells['rho'] * snap.info['unit_density'].express(C.g_cc)
    # for i in range(n_dust):
    #     output[f'densd{(i + 1):02d}'] = cells[f'rhod{(i + 1):02d}'] * snap.info['unit_density'].express(C.g_cc)

    # cell lengths
    unit_l = snap.info['unit_length'].express(C.m)
    output['dx'] = cells.get_sizes() * unit_l

    # max. number of cells
    numcells = len(output['dx'])
    
    # original cell positions (from 0 to 1) converted into uint length
    output['x'] = cells.points[:, 0] * unit_l
    output['y'] = cells.points[:, 1] * unit_l
    output['z'] = cells.points[:, 2] * unit_l
    print('\nx min max', output['x'].min(), output['x'].max())

    # level of each cell
    output['level'] = np.log2(unit_l / output['dx'])

    # gas temperature in K
    print('*************************************************')
    # ============================================
    # everything in cgs for the moment:
    if 'mu_gas' in snap.info:
        mu = snap.info['mu_gas'] # mean molecular mass in amu
    else:
        mu = 1.4
    print('WORKING WITH mu_gas = ', mu)
    # G        = 6.7e-08
    # kbol     = 1.38062e-16          # erg/degre
    # pc       = 3.08e+18             # cm
    # mp       = mu * 1.660531e-24    # n gramme
    G        = co.G.cgs.value
    kbol     = co.k_B.cgs.value
    pc       = co.pc.cgs.value
    mp       = mu * co.u.cgs.value
    scale_n  = 1.
    scale_l  = pc
    scale_d  = scale_n * mp
    scale_t  = 1.0 / np.sqrt(G * scale_d)
    scale_v  = scale_l / scale_t    
    scale_T2 = mp / kbol * scale_v**2
    #============================================
    unit = snap.info['unit_temperature'].express(C.K)
    print('unit', unit)
    print('scale_T2', scale_T2)
    # X = 0.76 # hydrogen mass fraction
    output['Tgas'] = cells['P'] / cells['rho'] * scale_T2
    output['Tdust'] = cells['P'] / cells['rho'] * scale_T2
    print('Min Max Mean T', output['Tgas'].min(), output['Tgas'].max(), output['Tgas'].mean())
    print('*************************************************')

    # velocity in cm/s
    output['vel'] = cells['vel'] * snap.info['unit_velocity'].express(C.cm / C.s)
    
    # B-field in G
    unit_b = snap.info['unit_mag'].express(C.Gauss)
    unit_b = unit_b * 8 * np.pi	
    output['B'] = 0.5 * (cells['Br'] + cells['Bl']) * unit_b

    return output, numcells, unit_l, n_dust

def convert_ramses2polaris(
        datapath: str,
        num: int,
        outpath: str,
        # star_positions: list = None,
        # hole_radius_au = 0,
        has_dust_in_sim = True
    ):

    global nr_of_cells
    global cell_counter

    o_num    = str(num).zfill(5)
    out_name = 'output_' + o_num

    input_file  = datapath + out_name
    output_file = outpath + 'ramses_grid_' + o_num + '.dat'
    
    # print ('input:', input_file)
    # print ('output:', output_file)
    # print ('Loading RAMSES data from: \n', input_file, '\n\n')
    
    # read RAMSES data
    data, nr_of_cells, max_length, n_dust =\
    loadRamsesData(input_file, has_dust_in_sim = has_dust_in_sim)

    # data IDs for the grid header
    data_ids = [4, 5, 6, 28, 3] # Bx, By, Bz, gas mass density, gas temperature
    data_ids += [29] * n_dust   # dust mass density

    # transpose center of cube
    x_min = -0.5 * max_length
    y_min = -0.5 * max_length
    z_min = -0.5 * max_length
    
    max_level = max(data['level'])
    min_level = min(data['level'])
    
    print ('\n\n\n')
    print ('Octree parameter:')
    print ('    Level        (min, max) :', int(min_level), ',' , int(max_level))
    print ('    Nr. of cells (data, max):', nr_of_cells, ',',  8**max_level)
    print ('    Length       (min, max) :', max_length / (2**max_level) , ',', max_length, '\n')
    
    # init. octree
    tree = OcTree(x_min, y_min, z_min, max_length)
            
    # fill octree
    for i in range(0, nr_of_cells):
        
        # create single cell
        level = data['level'][i]
        
        c_x = data['x'][i] - 0.5 * max_length
        c_y = data['y'][i] - 0.5 * max_length
        c_z = data['z'][i] - 0.5 * max_length

        Tgas = data['Tgas'][i]

        mag_x, mag_y, mag_z = data['B'][i]
        vel_x, vel_y, vel_z = data['vel'][i]
      
        # # If a star position is given, create a hole around it.
        # if star_positions:
        #     hole_radius_m = hole_radius_au * co.au.value
        #     for star_pos in star_positions:
        #         # Calculate distance from the cell to the current star.
        #         dstar = np.sqrt((c_x - star_pos[0])**2 + (c_y - star_pos[1])**2 + (c_z - star_pos[2])**2)
        #         # If the cell is within the hole radius, set densities to zero.
        #         if (dstar <= hole_radius_m):
        #             print("\npos", c_x, c_y, c_z)
        #             data['densgas'][i] = 0
        #             for j in range(n_dust): data[f'densd{(j + 1):02d}'][i] = 0
        #             print("dust mass density used to be", data['densd01'][i])

        densgas = data['densgas'][i]
        densdus = [data[f'densd{(j + 1):02d}'][i] for j in range(n_dust)]

        cell = cell_oct(c_x, c_y, c_z, 0, level)
        cell.data = [mag_x, mag_y, mag_z, densgas, Tgas] + densdus
        
        if i % 10000 == 0:
            sys.stdout.write('Constructing octree: ' + str(100.0 * i / nr_of_cells) + ' %    \r')
            sys.stdout.flush()
        
        # insert single cell into octree
        cell_root = tree.root
        tree.insertInTree(cell_root, cell, 0)
    
    sys.stdout.write(CLR_LINE)
    print ("Constructing octree:    done   ")
    
    # check octree integrity
    check = tree.checkOcTree(cell_root)
    sys.stdout.write(CLR_LINE)
        
    if check == False:
        print("ERROR: Octree integrity is inconsistent!   \n\n")
        exit()
    else:
        print ("Octree structure   :    OK      ")
    
    # write octree file header
    data_len = len(data_ids)
    file = open(output_file, 'wb')
        
    file.write(struct.pack('H', grid_id))
    file.write(struct.pack('H', data_len))

    for d_ids in data_ids:
        file.write(struct.pack('H', d_ids))

    file.write(struct.pack('d', max_length))
    
    # write octree
    cell_counter = 0.0
    tree.writeOcTree(file, tree.root)
    sys.stdout.write(CLR_LINE)

    print ("Writing octree     :    done   \n")
    print ("Octree successfully created")

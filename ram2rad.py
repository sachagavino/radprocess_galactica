# ram2rad.py

# Ilseung Han
#   27.07.2025
#   31.07.2025
#   13.08.2025
#   23.08.2025
#   01.09.2025
#   11.09.2025
#   30.10.2025
#   27.11.2025 (RAMES temperature)
#   14.12.2025
#   09.02.2026
#   11.03.2026

import os
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
grid_id = 20         # grid ID (20 = octree)  
# gamma   = 1.4        # Ramses default
# mH      = 1.66e-24   # Hydrogen mass in g
# kB      = 1.38e-16   # Boltzmann constant in cgs

CLR_LINE =   "                                                      \r"

max_level    = 0
cell_counter = 0
nr_of_cells  = 0

class cell_oct:
    """written by Stefan Reiss"""
    def __init__(self, _x_min, _y_min, _z_min, _length, _level):
        self.x_min = _x_min
        self.y_min = _y_min
        self.z_min = _z_min
        
        self.length = _length
        self.level = _level
    
        self.isleaf = 0
        self.data = []
        self.branches = []  
                  
class OcTree:
    """written by Stefan Reiss"""
    def __init__(self, _x_min, _y_min, _z_min, _length):
        self.root = cell_oct(_x_min, _y_min, _z_min, _length, 0)

    def initCellBoundaries(self, cell, _level):
        x_min = cell.x_min
        y_min = cell.y_min
        z_min = cell.z_min
        l = 0.5 * cell.length

        level = _level

        cell.isleaf = 0
        cell.data = []
        cell.branches = [None, None, None, None, None, None, None, None]
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
            cell_pos.data = cell.data
            cell_pos.isleaf = 1
                        
        else:    
            if len(cell_pos.branches) == 0:
                self.initCellBoundaries(cell_pos, _level + 1)

            x_mid = cell_pos.x_min + 0.5 * cell_pos.length
            y_mid = cell_pos.y_min + 0.5 * cell_pos.length
            z_mid = cell_pos.z_min + 0.5 * cell_pos.length
            
            new_cell_pos = cell_pos

            if(z_pos < z_mid): # z 0 1 2 3

                if(y_pos < y_mid): # y 0 1

                    if(x_pos < x_mid): # x 0
                        new_cell_pos = cell_pos.branches[0]
                    else: # x 1
                        new_cell_pos = cell_pos.branches[1]

                else: # y 2 3

                    if(x_pos < x_mid): # x 2
                        new_cell_pos = cell_pos.branches[2]
                    else: # x 3
                        new_cell_pos = cell_pos.branches[3]

            else: # z 4 5 6 7

                if(y_pos < y_mid): # y 4 5

                    if(x_pos < x_mid): # x 4
                        new_cell_pos = cell_pos.branches[4]
                    else: # x 5
                        new_cell_pos = cell_pos.branches[5]

                else: # y 6 7

                    if(x_pos < x_mid): # x 6
                        new_cell_pos = cell_pos.branches[6]
                    else: # x 7
                        new_cell_pos = cell_pos.branches[7]

            self.insertInTree(new_cell_pos, cell, _level + 1)

    def writeOcTree(self, file, cell):
        global cell_counter
        global nr_of_cells

        file.write(struct.pack('H', cell.isleaf))
        file.write(struct.pack('H', cell.level))   

        if cell.isleaf == 1:    
            data_len = len(cell.data)
            
            if cell_counter % 10000 == 0:
                sys.stdout.write("-> Writing octree grid file : "
                                 + str(100.0 * cell_counter / nr_of_cells)
                                 + " %     \r")
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
                sys.stdout.write("-> Checking octree integrity : "
                                 + str(100.0 * cell_counter / nr_of_cells)
                                 + " %     \r")
                sys.stdout.flush()
                
            cell_counter += 1    
            
        else:
            length = len(cell.branches)
            
            if length == 0:
                return False
            
            for i in range(8):
                self.checkOcTree(cell.branches[i])                
                
        return True

    # RADMC-3D                
    def writeOcTree_radmc(self, cell, grid, value):
        global cell_counter
        global nr_of_cells

        if cell.isleaf == 1:    
            data_len = len(cell.data)
            
            if cell_counter % 10000 == 0:
                sys.stdout.write("-> Writing octree grid file : "
                                 + str(100.0 * cell_counter / nr_of_cells)
                                 + " %     \r")
                sys.stdout.flush()
                
            cell_counter += 1
                
            value.append(cell.data)
            grid.append(0)

        else:
            grid.append(1)
            
            for i in range(8):
                self.writeOcTree_radmc(cell.branches[i], grid, value)
        
def loadRamsesData(filename, has_dust_in_sim = True):
    """written by Valeska Valdivia"""
    # Filename will be the format my/folder/here/output_?????
    # where ????? is from 00001 to 99999
    # print(filename)
    
    # Split the filename into folder and number (pymses needs this)
    outstr = 'output_'
    outloc = filename.rfind(outstr)
    folder = filename[:outloc]
    numloc = outloc + len(outstr)
    
    # Note: 5 characters in output number
    num = int(filename[numloc:numloc+5])

    # Create the pymses RamsesOutput object 
    # print("\n\nfolder: ", folder , num, "\n\n")
    snap = pymses.RamsesOutput(folder,num)

    # Spatial information
    output = {}

    if has_dust_in_sim:
        print("\n--- Loading data for a simulation WITH dust fields ---")
        # Load dust-to-fluid mass ratios for all available dust species
        # print(snap.info)
        n_dust = snap.info['ndust']
        dust_ratios = [f'dust_ratio_{i + 1}' for i in range(n_dust)]

        # Create a flat structure with the snapshot's cell data in it
        # snap.amr_fields()
        amr = snap.amr_source(['rho', 'P', 'vel', 'Br', 'Bl'] + dust_ratios)
        # print(amr)
        cell_source = CellsToPoints(amr)
        cells = cell_source.flatten()
        # print("The B fields", cells['Br'], cells['Bl'], "end")
   
        # Calculate densities based on dust ratios.
        # Mass density of fluid (in cgs)
        output['denstot'] = cells['rho'] * snap.info['unit_density'].express(C.g_cc)
        # print("dens min, max", output['denstot'].min(), output['denstot'].max())
   
        # Mass density of dust and gas (in cgs)
        output['epsilon'] = sum(cells[dust_ratio] for dust_ratio in dust_ratios)
        output['densdto'] = output['denstot'] * output['epsilon']
        for i in range(n_dust):
            output[f'densd{(i + 1):02d}'] = output['denstot'] * cells[f'dust_ratio_{i + 1}']
            # output[f'densd{(i + 1):02d}'] = output['densdto'] / n_dust # for test
        output['densgas'] = output['denstot'] - output['densdto']

    else:
        print("\n--- Loading data for a simulation WITHOUT dust fields ---")
        # Assume a single, virtual dust species.
        n_dust = 1
        amr = snap.amr_source(['rho', 'P', 'vel', 'Br', 'Bl'])
        cells = CellsToPoints(amr).flatten()

        # Calculate gas density directly and dust density as 1% of gas density.
        output['densgas'] = cells['rho'] * snap.info['unit_density'].express(C.g_cc)
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

    # Cell lengths
    unit_l = snap.info['unit_length'].express(C.m)
    output['dx'] = cells.get_sizes() * unit_l

    # Max. number of cells
    numcells = len(output['dx'])
    
    # Original cell positions (from 0 to 1) converted into uint length
    output['x'] = cells.points[:, 0] * unit_l
    output['y'] = cells.points[:, 1] * unit_l
    output['z'] = cells.points[:, 2] * unit_l
    print("\nx min max", output['x'].min(), output['x'].max())

    # Level of each cell
    output['level'] = np.log2(unit_l / output['dx'])

    # Gas temperature in K
    print("*************************************************")
    # ============================================
    # Everything in cgs for the moment:
    if 'mu_gas' in snap.info:
        mu = snap.info['mu_gas'] # mean molecular weight in amu
    else:
        mu = 1.4
    print("WORKING WITH mu_gas = ", mu)
    # G        = 6.7e-8
    # kbol     = 1.38062e-16          # erg/degre
    # pc       = 3.08e18              # cm
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
    scale_T2 = mp/kbol * scale_v**2
    # ============================================
    unit = snap.info['unit_temperature'].express(C.K)
    print("unit", unit)
    print("scale_T2", scale_T2)
    # X = 0.76 # Hydrogen mass fraction
    output['Tgas'] = cells['P'] / cells['rho'] * scale_T2
    output['Tdust'] = cells['P'] / cells['rho'] * scale_T2
    print("Min Max Mean T", output['Tgas'].min(), output['Tgas'].max(), output['Tgas'].mean())
    print("*************************************************")

    for i in range(n_dust):
        output[f'tempd{(i + 1):02d}'] = output['Tdust']

    # Velocity in m/s
    output['vel'] = cells['vel'] * snap.info['unit_velocity'].express(C.m / C.s)
    
    # B-field in G
    unit_b      = snap.info['unit_mag'].express(C.T)
    unit_b      = unit_b * 8 * np.pi
    output['B'] = 0.5 * (cells['Br'] + cells['Bl']) * unit_b

    return output, numcells, unit_l, n_dust

def convert_ramses2radmc3d(
        datapath: str,
        num: int,
        outpath: str,
        # star_positions: list = None,
        # hole_radius_au = 0,
        has_dust_in_sim = True
    ):

    """
    Function to call loadRamsesData() to load RAMSES outputs
    and write 'dust_density.inp' and 'amr_grid.inp' for radmc3d

    INPUTS:
        datapath (string)   : Path to the simulations data (that contains different output folders)
        num (int)           : Simulation output number
        dust_to_gas (float) : Dust to gas mass ratio
        outpath (string)    : Directory to write output files to

    OUTPUT:
        'dust_density.inp' and 'amr_grid.inp' files in outpath
    """

    global nr_of_cells
    global cell_counter
   
    # print(pymses.__file__)

    # print ("====================================================================================================")
    o_num      = str(num).zfill(5)
    out_name   = 'output_' + o_num
    input_file = datapath + out_name

    # print ("Ramses input:    ", input_file)
    # print ("Output for RADMC: amr_grid.inp, dust_density.inp")
    # print ("")

    # print ("Loading RAMSES data from: \n", input_file)

    # Read RAMSES data
    data, nr_of_cells, max_length, n_dust =\
    loadRamsesData(input_file, has_dust_in_sim = has_dust_in_sim)
    L_cm = pymses.RamsesOutput(datapath, int(num)).info["unit_length"].express(C.cm) # box length in cm

    # Transpose center of cube
    x_min = -0.5 * max_length
    y_min = -0.5 * max_length
    z_min = -0.5 * max_length
    
    max_level = max(data['level'])
    min_level = min(data['level'])
    
    print("\n")
    print("Octree parameter:")
    print("    Level        (min,max)  : ", int(min_level), ",", int(max_level))
    print("    Nr. of cells (data, max): ", nr_of_cells, ",", 8**max_level)
    print("    Length       (min,max)  : ", max_length / (2**max_level), ",", max_length, "\n")
    
    # Init. octree
    tree = OcTree(x_min, y_min, z_min, max_length)

    # Fill octree (for grid and density)
    for i in range(0, nr_of_cells):
        # Create single cell
        level = data['level'][i]
        
        c_x = data['x'][i] - 0.5 * max_length
        c_y = data['y'][i] - 0.5 * max_length
        c_z = data['z'][i] - 0.5 * max_length

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
        #             # data['densgas'][i] = 0
        #             for j in range(n_dust): data[f'densd{(j + 1):02d}'][i] = 0
        #             print("dust mass density used to be", data['densd01'][i])

        densdus = [data[f'densd{(j + 1):02d}'][i] for j in range(n_dust)]

        cell = cell_oct(c_x, c_y, c_z, 0, level)
        cell.data = densdus

        # Insert single cell into octree
        cell_root = tree.root
        tree.insertInTree(cell_root, cell, 0)

        if i % 10000 == 0:
            sys.stdout.write("Constructing octree: " + str(100.0 * i / nr_of_cells) + " %    \r")
            sys.stdout.flush()

    sys.stdout.write(CLR_LINE)
    print ("\nConstructing octree:    done   ")

    # Check octree integrity
    print ("Calling tree.checkOcTree(cell_root), nr_of_cells = ", nr_of_cells)
    check = tree.checkOcTree(cell_root)
    # print ("Tree OK ;)")
    sys.stdout.write(CLR_LINE)
        
    if check == False:
        print ("ERROR: Octree integrity is inconsistent!   \n\n")
        exit ()
    else:
        print ("Octree structure   :    OK      ") 
            
    print("Writing octree     :    done   \n")
    
    print("Octree successfully created\n")

    # Write octree
    cell_counter = 0.0    
    grid         = []   # Vector containing cell values (0/1) in RADMC convention
    density      = []   # Vector containing dust density in RADMC convention
    tree.writeOcTree_radmc(tree.root, grid, density)
    densityarray = np.array(density)
    sys.stdout.write(CLR_LINE)

    # Create output directory if not existed
    os.makedirs(os.path.dirname(outpath), exist_ok = True)
    print("Writing the amr_grid.inp file for RADMC-3D...\n")
    with open(outpath + 'amr_grid.inp','w+') as f:
        f.write("1\n")                                                  # iformat (typically 1 at present)
        f.write("1\n")                                                  # grid style (1: Oct-tree)
        f.write("1\n")                                                  # coordinates (cartesian if < 100)
        f.write("0\n")                                                  # gridinfo (0 recommended)
        f.write("1\t1\t1\n")                                            # incl_x, incl_y, incl_z
        f.write("1\t1\t1\n")                                            # nx, ny, nz
        f.write("%d\t%d\t%d\n" % (max_level, nr_of_cells, len(grid)))   # levelmax, nleafsmax, nbranchmax
        f.write("%e\t%e\n" % (-L_cm / 2, L_cm / 2))
        f.write("%e\t%e\n" % (-L_cm / 2, L_cm / 2))
        f.write("%e\t%e\n" % (-L_cm / 2, L_cm / 2))

        for i in range(len(grid)):
            f.write("%d\n" % grid[i])
    
    print("Writing the dust_density.inp file for RADMC-3D...")
    with open(outpath + 'dust_density.inp','w+') as f:
        f.write("1\n")
        f.write("%d\n" % nr_of_cells)
        f.write("%d" % densityarray.shape[1])
       
        for i in range(densityarray.shape[1]):
            for j in range(densityarray.shape[0]):
                f.write("\n%e" % densityarray[j, i])
        # return density, densityarray

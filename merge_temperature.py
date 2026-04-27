# merge_temperature.py

# Ilseung Han
#   21.08.2025
#   23.08.2025
#   01.09.2025
#   04.01.2026

import os
import yaml
import struct
import numpy as np

from utils import check_simulation_has_dust, get_dust_species_count

def main():
    """
    Main function to merge RADMC-3D temperature data into a POLARIS grid file.
    This script renames the original POLARIS grid to 'grid_temp.polaris.dat'
    and saves the new, combined grid as 'grid_temp.radmc3d.dat'.
    """
    config_file_path = 'config.yaml'
    try:
        with open(config_file_path, 'r') as f:
            config = yaml.safe_load(f)
            print(f"\nLoaded configuration from: '{config_file_path}'")

        polaris_run_dir = config['polaris_opacity_run']['output_path']
        radmc3d_run_dir = config['radmc3d_run']['run_directory']

        # Define the paths for the files within the POLARIS run directory.
        original_grid = os.path.join(polaris_run_dir, 'grid_temp.dat')
        renamed_polaris_grid = os.path.join(polaris_run_dir, 'grid_temp.polaris.dat')
        final_radmc3d_grid = os.path.join(polaris_run_dir, 'grid_temp.radmc3d.dat')

        # Define the path to the temperature data from the RADMC-3D run.
        radmc3d_temp_file = os.path.join(radmc3d_run_dir, 'dust_temperature.bdat')

        # Rename the original POLARIS grid file.
        if os.path.exists(original_grid):
            # Check if the renamed file already exists to prevent errors on re-runs.
            if os.path.exists(renamed_polaris_grid):
                os.remove(renamed_polaris_grid)
            os.rename(original_grid, renamed_polaris_grid)
            print(f"\nRenamed '{original_grid}' to '{renamed_polaris_grid}'")
        elif not os.path.exists(renamed_polaris_grid):
            raise FileNotFoundError(f"Neither '{original_grid}' nor '{renamed_polaris_grid}' were found.")
        else:
            print(f"\nOriginal file '{renamed_polaris_grid}' already exists. Using it directly.")

        ramses_dir = config['ramses_output_dir']
        output_num = int(os.path.basename(os.path.normpath(ramses_dir)).split('_')[-1])

        info_file = os.path.join(ramses_dir, f"info_{output_num:05d}.txt")
        has_dust = check_simulation_has_dust(info_file)

        if has_dust:
            n_dust = get_dust_species_count(info_file)
        else:
            n_dust = 1
            # print("  - Assuming a single dust species derived from gas density.")
            print("\nAssuming a single dust species derived from gas density.")

        # Load the RADMC-3D temperature data.
        with open(radmc3d_temp_file, 'rb') as f:
            f.seek(4 * 8) # Skip the 4-header entries (each 8 bytes).
            temperatures = np.fromfile(f, dtype=np.float32)

        # Reshape the 1D temperature array into a 2D array of (n_cells, n_dust_species).
        temperatures_reshaped = temperatures.reshape(-1, n_dust, order = 'F')
        num_leaf_cells = temperatures_reshaped.shape[0]
        print(f"\nLoaded and reshaped temperatures for {num_leaf_cells} cells and {n_dust} dust species.")

        print(f"\nMerging data and writing new grid file to '{final_radmc3d_grid}'")

        leaf_index = 0
        # Copy header information directly from the original POLARIS grid.
        with open(renamed_polaris_grid, 'rb') as in_f,\
             open(final_radmc3d_grid, 'wb') as out_f:
            grid_id = in_f.read(2)
            out_f.write(grid_id)

            n_param_buf = in_f.read(2)
            out_f.write(n_param_buf)
            n_param = struct.unpack('H', n_param_buf)[0]

            # Read all parameter IDs into a list.
            param_ids_buf = in_f.read(2 * n_param)
            param_ids = list(struct.unpack('H' * n_param, param_ids_buf))
            out_f.write(param_ids_buf) # Write the original IDs to the new file.

            # Find all indices for dust temperature (ID = 2).
            temp_indices = [i for i, pid in enumerate(param_ids) if pid == 2]

            if len(temp_indices) != n_dust:
                raise ValueError(f"Mismatch: Found {len(temp_indices)} temperature parameters (ID = 2),\
                                   but expected {n_dust} dust species.")
            print(f"  - Found {len(temp_indices)} dust temperature parameters at indices: {temp_indices}")

            grid_size = in_f.read(8)
            out_f.write(grid_size)

            # Loop through the grid structure until all leaf cells are processed.
            while leaf_index < num_leaf_cells:
                buf = in_f.read(2)
                if not buf: break
                out_f.write(buf)
                is_leaf = struct.unpack('H', buf)[0]

                buf = in_f.read(2)
                if not buf: break
                out_f.write(buf)

                if is_leaf == 1:
                    dust_temp_counter = 0
                    for j in range(n_param):
                        val_buf = in_f.read(4)
                        # If the current parameter index is one of the temperature indices...
                        if j in temp_indices:
                            # ...replace it with the corresponding temperature from RADMC-3D.
                            temp_val = temperatures_reshaped[leaf_index, dust_temp_counter]
                            out_f.write(struct.pack('f', temp_val))
                            dust_temp_counter += 1
                        else:
                            # Otherwise, just copy the original value.
                            out_f.write(val_buf)
                    leaf_index += 1
        
        print(f"Processed {leaf_index} leaf cells.")
        print(f"\nNew POLARIS grid with RADMC-3D temperatures saved successfully!")

    except Exception as e:
        print(f"\nERROR in pipeline: {e}")

if __name__ == '__main__':
    main()

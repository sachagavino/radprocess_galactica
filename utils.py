# utils.py

# Ilseung Han
#   19.08.2025
#   24.08.2025
#   31.08.2025
#   01.09.2025
#   02.09.2025
#   03.09.2025
#   04.09.2025
#   10.09.2025
#   11.09.2025
#   14.09.2025
#   26.09.2025
#   13.10.2025
#   30.10.2025
#   31.10.2025
#   14.12.2025
#   19.01.2026

import os
import numpy as np
import astropy.units as u
import astropy.constants as co

def check_simulation_has_dust(info_file_path: str) -> bool:
    """
    Checks if a RAMSES simulation contains dust
    by looking for the 'ndust' parameter in the info file.

    Args:
        info_file_path (str): The path to the info_xxxxx.txt file.
    Returns:
        bool: True if 'ndust' is found, False otherwise.
    """
    print(f"\nChecking 'ndust' for dust fields in: '{os.path.basename(info_file_path)}'")
    try:
        with open(info_file_path, 'r') as f:
            for line in f:
                # Check if the line starts with 'ndust' after stripping whitespace.
                if line.strip().startswith('ndust'):
                    print("  - Found 'ndust' parameter. Assuming simulation has dust.")
                    return True
            # If the loop finishes without finding 'ndust', it doesn't have dust.
        print("  - Not found 'ndust' parameter. Assuming simulation is gas-only.")
        return False
    
    except FileNotFoundError:
        # Re-raise the error to be caught by the main script's try/except block.
        # This informs the calling script that the file is missing and the pipeline should halt.
        raise FileNotFoundError(f"Info file not found at: {info_file_path}")

def derive_stars_properties(ramses_output_dir: str) -> list[dict]:
    """
    Derives comprehensive stellar properties (L, Teff, R) for all stars
    based on their mass from the sink file and physical assumptions.
    """
    all_stars_basic_info = get_stars_properties(ramses_output_dir)
    if not all_stars_basic_info:
        return []

    derived_stars_list = []
    for star in all_stars_basic_info:
        # if star['lint'] > 0:
        if star['lint'] != 0:
            ltot = star['lacc'] + star['lint']
            ltot = ltot * u.L_sun
            teff = star['teff'] * u.K
            radi = np.sqrt(ltot / (4 * np.pi * co.sigma_sb * teff**4))  # Tung et al. (2025)              
            radi = radi.to(u.R_sun)
            teff = teff.to(u.K)

            ltot = ltot.value
            teff = teff.value
            radi = radi.value

            derived_star_props = {
                'id': star['id'],
                'mass_msun': star['mass'],
                'pos_m': star['pos'],
                'radius_rsun': radi,
                'luminosity_lsun': ltot,
                'temperature_k': teff
            }
            derived_stars_list.append(derived_star_props)

    return derived_stars_list

def get_dust_species_count(info_file_path: str) -> int:
    """
    Reads the number of dust species from the RAMSES info file.

    Args:
        info_file_path (str): The path to the info_xxxxx.txt file.

    Returns:
        int: The number of dust species found.
    """
    print(f"\nReading info file to find 'ndust': {os.path.basename(info_file_path)}")
    try:
        with open(info_file_path, 'r') as f:
            for line in f:
                if 'ndust' in line.strip():
                    n_dust = int(line.split('=')[1])
                    print(f"  - Found {n_dust} dust species.")
                    return n_dust
        # If the loop finishes without finding 'ndust', raise an error.
        raise ValueError(f"'ndust' key not found in {info_file_path}")

    except FileNotFoundError:
        # Re-raise the error to be caught by the main function's try/except block.
        raise FileNotFoundError(f"Info file not found at: {info_file_path}")

def get_sink_format(ramses_output_dir: str, output_num_str: str) -> dict:
    """
    Automatically determines the column format of the sink file
    by checking for a sink.info file or a header in the sink.csv file.
    """
    print("\nDetermining sink file format...")
    sink_info_file = os.path.join(ramses_output_dir, f'sink_{output_num_str}.info')
    sink_csv_file = os.path.join(ramses_output_dir, f'sink_{output_num_str}.csv')

    # Priority 1: Check for sink.info file.
    if os.path.exists(sink_info_file):
        with open(sink_info_file, 'r') as f:
            for line in f:
                if 'M[Msol]' in line and 'x' in line:
                    headers = line.strip().split()
                    try:
                        return {
                            'id_col': headers.index('Id'),
                            'mass_col': headers.index('M[Msol]'),
                            'x_col': headers.index('x'),
                            'y_col': headers.index('y'),
                            'z_col': headers.index('z'),
                            'lacc_col': headers.index('acc_lum[Lsol]'),
                            'lint_col': headers.index('int_lum[Lsol]'),
                            'teff_col': headers.index('Teff[K]')
                        }
                    except ValueError:
                        continue

    # Priority 2: Check for a header in sink.csv
    if os.path.exists(sink_csv_file):
        with open(sink_csv_file, 'r') as f:
            first_line = f.readline().strip().lstrip('#')
            if any(c.isalpha() for c in first_line):
                headers = [h.strip() for h in first_line.split(',')]
                try:
                    return {
                        'id_col': headers.index('id'),
                        'mass_col': headers.index('msink'),
                        'x_col': headers.index('x'),
                        'y_col': headers.index('y'),
                        'z_col': headers.index('z'),
                    }
                except ValueError:
                    pass

    raise ValueError(f"Could not automatically determine sink file format for output {output_num_str}.")

def get_star_positions(ramses_output_dir: str) -> list[list[float]]:
    """
    A convenience function to get a list of position for all stars.
    """
    all_stars = get_stars_properties(ramses_output_dir)
    positions = [star['pos'] for star in all_stars]
    
    return positions

def get_stars_properties(ramses_output_dir: str) -> list[dict]:
    """
    Reads RAMSES output files to get the properties of all stars (sink particles).
    """
    # print("\nReading properties for all stars from:")
    # print(f"  '{ramses_output_dir}/sink_{output_num_str}.csv'")
    
    try:
        output_num_str = os.path.basename(os.path.normpath(ramses_output_dir)).split('_')[-1]
        sink_format = get_sink_format(ramses_output_dir, output_num_str)

        info_file = os.path.join(ramses_output_dir, f'info_{output_num_str}.txt')
        sink_file = os.path.join(ramses_output_dir, f'sink_{output_num_str}.csv')

        print("\nReading properties for all stars from:")
        print(f"  '{ramses_output_dir}/sink_{output_num_str}.csv'")

        with open(info_file, 'r') as f:
            box_size_pc = 0.0
            length_unit_cm = 0.0
            for line in f:
                if 'boxlen' in line:
                    box_size_pc = float(line.split('=')[1])
                    # print(f"  - Found box size in info file: {box_size_pc} pc")
                    print(f"  - Found box size in 'info_{output_num_str}.txt': {box_size_pc} pc")
                elif 'unit_l' in line:
                    length_unit_cm = float(line.split('=')[1])
                    # print(f"  - Found length unit in info file: {length_unit_cm} cm")
                    print(f"  - Found length unit in 'info_{output_num_str}.txt': {length_unit_cm} cm")
        length_unit_m = length_unit_cm * 1e-2

        sink_data = np.loadtxt(sink_file, delimiter = ',', ndmin = 2, comments = ' #')

        stars_list = []
        for star_row in sink_data:
            pos_pc = [
                star_row[sink_format['x_col']],
                star_row[sink_format['y_col']],
                star_row[sink_format['z_col']]
            ]
            # print(f"\nFound star position in sink file: {pos_pc} pc")
            print(f"\nFound star position in 'sink_{output_num_str}.csv': {pos_pc} pc")
            pos_m = [(p - box_size_pc / 2.0) * length_unit_m for p in pos_pc]
            print(f"Final star position: {pos_m} m")
        
            stars_list.append({
                'id': int(star_row[sink_format['id_col']]),
                'mass': star_row[sink_format['mass_col']],
                'pos': pos_m,
                'lacc': star_row[sink_format['lacc_col']],
                'lint': star_row[sink_format['lint_col']],
                'teff': star_row[sink_format['teff_col']]
            })

        print(f"\nFound and processed {len(stars_list)} stars.")
        return stars_list

    except Exception as e:
        print(f"ERROR: Failed to process star properties: {e}")
        # return []
        raise

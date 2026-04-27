# create_polaris_grid.py

# Ilseung Han
#   17.08.2025
#   18.08.2025
#   20.08.2025
#   23.08.2025
#   24.08.2025
#   31.08.2025
#   30.10.2025
#   31.10.2025
#   19.01.2026
#   11.03.2026

import os
import yaml

from ram2pol import convert_ramses2polaris
from utils import check_simulation_has_dust, get_star_positions

def main():
    """
    Main function for the POLARIS grid creation step.
    """
    config_file_path = 'config.yaml'
    try:
        with open(config_file_path, 'r') as f:
            config = yaml.safe_load(f)
            print(f"\nLoaded configuration from: '{config_file_path}'")

        ramses_dir = config['ramses_output_dir']
        output_num = int(os.path.basename(os.path.normpath(ramses_dir)).split('_')[-1])
        info_file = os.path.join(ramses_dir, f'info_{output_num:05d}.txt')
        has_dust = check_simulation_has_dust(info_file)

        # star_pos = get_star_positions(ramses_dir)
        # if not star_pos:
        #     raise ValueError("\nCould not calculate star position.")

        grid_config = config['grid_conversion']
        polaris_config = grid_config['polaris']
        output_path = polaris_config['output_path']
        # hole_radius_au = grid_config.get('hole_radius_au', 0)

        print("\nStarting RAMSES to POLARIS grid conversion...")
        os.makedirs(output_path, exist_ok = True)

        convert_ramses2polaris(
            datapath = os.path.dirname(ramses_dir),
            num = output_num,
            outpath = os.path.join(output_path, ''),
            # star_positions = star_pos,
            # hole_radius_au = hole_radius_au,
            has_dust_in_sim = has_dust
        )

        print("\nPOLARIS grid file created successfully!")

    except FileNotFoundError:
        print(f"\nERROR: Configuration file not found at '{config_file_path}'")
        raise
    except KeyError as e:
        print(f"\nERROR: Missing key in '{config_file_path}': {e}")
        raise
    except Exception as e:
        print(f"\nERROR: An unexpected error occurred in the pipeline: {e}")
        raise

if __name__ == '__main__':
    main()

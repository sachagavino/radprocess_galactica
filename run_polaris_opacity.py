# run_polaris_opacity.py

# Ilseung Han
#   18.08.2025
#   19.08.2025
#   23.08.2025
#   01.09.2025
#   02.09.2025
#   03.09.2025
#   30.10.2025
#   31.10.2025
#   18.01.2026
#   19.01.2026
#   09.02.2026
#   12.03.2026
#   13.03.2026

import os
import sys
import yaml
import shutil
import subprocess
import numpy as np

from utils import check_simulation_has_dust, derive_stars_properties, get_dust_species_count

def cleanup_previous_run(output_path: str):
    """
    Safely removes output files and directories from a previous POLARIS run.
    """
    print("\nCleaning up previous run outputs...")

    data_dir = os.path.join(output_path, 'data')
    plots_dir = os.path.join(output_path, 'plots')
    grid_temp_file = os.path.join(output_path, 'grid_temp.dat')

    # shutil.rmtree removes a directory and everything inside it.
    if os.path.exists(data_dir):
        shutil.rmtree(data_dir)
        print(f"  - Removed directory: '{data_dir}'")

    if os.path.exists(plots_dir):
        shutil.rmtree(plots_dir)
        print(f"  - Removed directory: '{plots_dir}'")

    # os.remove removes a single file.
    if os.path.exists(grid_temp_file):
        os.remove(grid_temp_file)
        print(f"  - Removed file: '{grid_temp_file}'")

def get_mf_grain_sizes(ramses_dir: str, output_num: int):
    filename = f'info_mf_{output_num:05d}.txt'
    filepath = os.path.join(ramses_dir, filename)

    if not os.path.exists(filepath):
        return []

    grain_sizes = []
    print(f"\nFound multi-fluid info file: '{filename}'")

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('mf_grain_size_'):
                    parts = line.split('=')
                    if len(parts) == 2:
                        try:
                            val = float(parts[1].strip())
                            grain_sizes.append(val)
                        except ValueError:
                            continue

        grain_sizes.sort()
        return grain_sizes

    except Exception as e:
        print(f"WARNING: Error reading {filename}: {e}")

def create_polaris_cmd_file(
        config: dict,
        all_stars_derived_props: list,
        # specific_grain_sizes: list,
        n_dust: int,
        output_num: int
    ):
    """
    Generates the POLARIS command file (.cmd) for the opacity calculation run.
    """
    polaris_config = config['polaris_opacity_run']
    output_path = polaris_config['output_path']
    cmd_filename = f'polaris_opacity_{output_num:05d}.cmd'
    cmd_filepath = os.path.join(output_path, cmd_filename)
    grid_input_file = os.path.join(polaris_config['grid_input_dir'], f'ramses_grid_{output_num:05d}.dat')

    os.makedirs(output_path, exist_ok = True)
    print(f"\nGenerating POLARIS command file at: '{cmd_filepath}'")

    with open(cmd_filepath, 'w') as f:
        # --- <common> block ---
        f.write('<common>\n')

        '''
        if specific_grain_sizes:
            print(f"\nUsing {len(specific_grain_sizes)} specific grain sizes from info file.")
            for i, size in enumerate(specific_grain_sizes):
                size_min = size * 0.9
                size_max = size * 1.1
                if i < 10:
                    f.write(f'\n\t<dust_component id = "{i}">  "{polaris_config["dust_nk_path"]}" "plaw" 1.0 0 '
                            f'{size_min:.4e} {size_max:.4e} 0.0')
                else:
                    f.write(f'\n\t<dust_component id = "{i}"> "{polaris_config["dust_nk_path"]}" "plaw" 1.0 0 '
                            f'{size_min:.4e} {size_max:.4e} 0.0')

        else:
            print(f"\nUsing settings from config.yaml.")
            size_bins = np.logspace(np.log10(polaris_config['dust_size_min']),
                                    np.log10(polaris_config['dust_size_max']), n_dust + 1)
            for i in range(n_dust):
                f.write(f'\n\t<dust_component id = "{i}"> "{polaris_config["dust_nk_path"]}" "plaw" 1.0 0 '
                        f'{size[i]:.2e} {size[i+1]:.2e} {polaris_config["dust_size_powerlaw"]}')
        '''
        print(f"\nUsing settings from config.yaml.")
        size = np.logspace(np.log10(polaris_config['dust_size_min']),
                           np.log10(polaris_config['dust_size_max']), n_dust + 1)
        for i in range(n_dust):
            # f.write(f'\n\t<dust_component id = "{i}"> "{polaris_config["dust_nk_path"]}" "plaw" 1.0 0 '
            #         f'{size[i]:.2e} {size[i+1]:.2e} {polaris_config["dust_size_powerlaw"]}')
            f.write(f'\n\t<dust_component id = "{i}"> "{polaris_config["dust_cs_path"][0]}" "plaw" 0.625 0 '
                    f'{size[i]:.2e} {size[i+1]:.2e} {polaris_config["dust_size_powerlaw"]}')
            f.write(f'\n\t<dust_component id = "{i}"> "{polaris_config["dust_cs_path"][1]}" "plaw" 0.375 0 '
                    f'{size[i]:.2e} {size[i+1]:.2e} {polaris_config["dust_size_powerlaw"]}')

        f.write(f'\n\n\t<nr_threads> {polaris_config["nr_threads"]}\n')
        f.write('\n</common>\n')

        # --- <task> block ---
        f.write('\n<task> 1\n')
        f.write('\n\t<cmd> CMD_TEMP\n')

        for star in all_stars_derived_props:
            xpos, ypos, zpos = star['pos_m']
            star_radius = star['radius_rsun']
            star_temp = star['temperature_k']
            # f.write(f'\n\t<source_star nr_photons = "1"> {xpos} {ypos} {zpos} {star_radius} {star_temp}')
            f.write(f'\n\t<source_star nr_photons = "1"> {xpos:17.10e} {ypos:17.10e} {zpos:17.10e} {star_radius:17.10e} {star_temp:17.10e}')
        # f.write('\n\t<source_isrf nr_photons = "0">\n')
        f.write('\n\t<source_isrf nr_photons = "1">  1\n')

        f.write(f'\n\t<path_grid> "{grid_input_file}"')
        f.write(f'\n\t<path_out>  "{os.path.join(output_path, "")}"\n')

        f.write(f'\n\t<mu> {polaris_config["mean_molecular_weight"]}\n')
        f.write(f'\n\t<mass_fraction> {polaris_config["mass_fraction"]}\n')
        f.write('\n</task>')

    return cmd_filepath

def main():
    """Main function to run the POLARIS opacity calculation."""
    config_file_path = 'config.yaml'
    try:
        with open(config_file_path, 'r') as f:
            config = yaml.safe_load(f)
            print(f"\nLoaded configuration from: '{config_file_path}'")

        ramses_dir = config['ramses_output_dir']
        output_num = int(os.path.basename(os.path.normpath(ramses_dir)).split('_')[-1])
        print(f"\nProcessing RAMSES output number: {output_num}")

        # all_stars = derive_stars_properties(ramses_dir)
        all_stars = derive_stars_properties('/data/pebbles/scratch/ilseung/pipeline/outputs/12-output_00013/')
        if not all_stars:
            raise ValueError("Could not calculate star position.")

        '''
        specific_grain_sizes = get_mf_grain_sizes(ramses_dir, output_num)
        n_dust = 0

        if specific_grain_sizes:
            n_dust = len(specific_grain_sizes)
        else:
            print(f"\n'info_mf' file not found.")
            info_file = os.path.join(ramses_dir, f"info_{output_num:05d}.txt")
            has_dust = check_simulation_has_dust(info_file)

            if has_dust:
                n_dust = get_dust_species_count(info_file)
            else:
                n_dust = 1
                print("\nAssuming a single dust species derived from gas density.")

        cmd_file = create_polaris_cmd_file(config, all_stars, specific_grain_sizes, n_dust, output_num)
        '''

        info_file = os.path.join(ramses_dir, f"info_{output_num:05d}.txt")                                   
        has_dust = check_simulation_has_dust(info_file)

        if has_dust:
            n_dust = get_dust_species_count(info_file)
        else:
            n_dust = 1
            print("\nAssuming a single dust species derived from gas density.")

        cmd_file = create_polaris_cmd_file(config, all_stars, n_dust, output_num)

        output_path = config['polaris_opacity_run']['output_path']
        cleanup_previous_run(output_path)
        
        log_filename = f'polaris_opacity_{output_num:05d}.log'
        log_filepath = os.path.join(output_path, log_filename)
        print(f"\nExecuting POLARIS... Log will be shown below and saved to: '{log_filepath}'\n")
        
        command = ['polaris', cmd_file]
        # command = ['/data/pebbles/software/polaris/github/bin/polaris', cmd_file]
        process = subprocess.Popen(
            command,
            stdout = subprocess.PIPE,   # Create a pipe for standard output.
            stderr = subprocess.STDOUT, # Redirect error messages to the same pipe.
            text = True,
            encoding = 'utf-8'
        )

        with open(log_filepath, 'w') as log_file:
            for line in process.stdout:
                sys.stdout.write(line)
                sys.stdout.flush()
                clean_line = line.rstrip('\r\n')
                log_file.write(clean_line + '\n')

        process.wait()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, command)

        print("POLARIS opacity run completed successfully!")

    except (ValueError, FileNotFoundError, KeyError) as e:
        print(f"\nERROR: Halting pipeline due to a configuration or file issue.")
        print(f"    Details: {e}")
        raise
    except subprocess.CalledProcessError as e:
        print(f"\nERROR: POLARIS exited with a non-zero status code: {e.returncode}")
        print(f"    Check the log file for details: {log_filepath}")
        raise
    except Exception as e:
        print(f"\nERROR: An unexpected error occurred in the pipeline: {e}")
        raise

if __name__ == '__main__':
    main()

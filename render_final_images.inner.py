# render_final_images.py

# Ilseung Han
#   22.08.2025
#   23.08.2025
#   01.09.2025
#   04.01.2026
#   14.03.2026

import os
import sys
import yaml
import shutil
import subprocess
import numpy as np
import astropy.units as u

from utils import check_simulation_has_dust, get_dust_species_count

def create_imaging_cmd_file(config: dict, output_num: int, n_dust: int, view_name: str, view_details: dict):
    """
    Generates a POLARIS command file for a specific viewing angle (e.g., 'xy', 'xz', 'yz').
    """
    opacity_config = config['polaris_opacity_run']
    image_config = config['final_image_rendering']
    output_base_path = image_config['output_base_path']
    # image_output_path = os.path.join(output_base_path, f'{view_name}/')
    image_output_path = os.path.join(output_base_path, f'inner/{view_name}/')
    
    # cmd_filename = f"image_{view_name}_{output_num:05d}.cmd"
    # cmd_filename = f'polaris_render_{view_name}_{output_num:05d}.cmd'
    base_filename = f'polaris_render_{view_name}_{output_num:05d}'
    # cmd_filename = f'{base_filename}.cmd'
    cmd_filename = f'{base_filename}.inner.cmd'
    cmd_filepath = os.path.join(output_base_path, cmd_filename)

    polaris_run_dir = config['polaris_opacity_run']['output_path']
    grid_input_file = os.path.join(polaris_run_dir, 'grid_temp.radmc3d.dat')

    print(f"\nGenerating POLARIS command file for '{view_name}' view: '{cmd_filepath}'")
    with open(cmd_filepath, 'w') as f:
        # --- <common> block ---
        f.write('<common>\n')

        size = np.logspace(np.log10(opacity_config['dust_size_min']),
                           np.log10(opacity_config['dust_size_max']), n_dust + 1)
        for i in range(n_dust):
            # f.write(f'\n\t<dust_component id = "{i}"> "{opacity_config["dust_nk_path"]}" "plaw" 1.0 0 '
            #         f'{size[i]:.2e} {size[i+1]:.2e} {opacity_config["dust_size_powerlaw"]}\n')
            f.write(f'\n\t<dust_component id = "{i}"> "{opacity_config["dust_cs_path"][0]}" "plaw" 0.625 0 '
                    f'{size[i]:.2e} {size[i+1]:.2e} {opacity_config["dust_size_powerlaw"]}')
            f.write(f'\n\t<dust_component id = "{i}"> "{opacity_config["dust_cs_path"][1]}" "plaw" 0.375 0 '
                    f'{size[i]:.2e} {size[i+1]:.2e} {opacity_config["dust_size_powerlaw"]}')
        # f.write(f'\n\n\t<nr_threads> {opacity_config["nr_threads"]}\n')
        f.write(f'\n\n\t<nr_threads> {image_config["nr_threads"]}\n')

        f.write('\n</common>\n')

        # --- <task> block ---
        f.write('\n<task> 1\n')

        f.write('\n\t<cmd> CMD_DUST_EMISSION\n')

        f.write(f'\n\t<path_grid> "{grid_input_file}"')
        # f.write(f'\n\t<path_out>  "{image_output_path}"\n')
        f.write(f'\n\t<path_out>  "{image_output_path}"\n')

        f.write(f'\n\t<mu> {opacity_config["mean_molecular_weight"]}\n')
        f.write(f'\n\t<mass_fraction> {opacity_config["mass_fraction"]}\n')
        
        f.write(f'\n\t<align> ALIG_PA\n')

        npix = image_config['npix']
        plane_id = view_details['plane_id']
        f.write(f'\n\t<write_inp_midplanes> {npix}')
        f.write(f'\n\t<write_3d_midplanes> {plane_id} {npix}\n')

        # f.write(f'\n\t<midplane_zoom> 1\n')
        f.write(f'\n\t<midplane_zoom> 10\n')

        axis1 = view_details['axis1']
        axis2 = view_details['axis2']
        f.write(f'\n\t<axis1> {axis1[0]} {axis1[1]} {axis1[2]}')
        f.write(f'\n\t<axis2> {axis2[0]} {axis2[1]} {axis2[2]}\n')

        # npix = image_config['npix']
        # wave_val = (image_config['wavelength_mm'] * u.mm).to(u.m).value
        dist_val = (image_config['distance_pc'] * u.pc).to(u.m).value
        # dist_val = image_config['distance_pc'] * constants['parsec_in_meters']
        boxl = 0.662166503408534E+02 * u.pc # orion/info_00013.txt
        size = boxl / 10
        size = size.to(u.m)
        size_val = size.value
        theta = view_details['theta']
        phi = view_details['phi']
        for wave_val in image_config['wavelength_mm']:
            wave_val = (wave_val * u.mm).to(u.m).value
            # f.write(f'\n\t<detector_dust nr_pixel = "{npix}"> {wave_val:e} {wave_val:e} 1 1 '
            #         f'{theta} {phi} {dist_val:e}')
            f.write(f'\n\t<detector_dust nr_pixel = "{npix}"> {wave_val:e} {wave_val:e} 1 1 '
                    f'{theta} {phi} {dist_val:e} {size_val:e} {size_val:e}')

        f.write('\n\n</task>')

    return cmd_filepath, image_output_path

def main():
    """Main function to render all final images with POLARIS."""
    config_file_path = 'config.yaml'
    try:
        with open(config_file_path, 'r') as f:
            config = yaml.safe_load(f)
            print(f"\nLoaded configuration from: '{config_file_path}'")

        ramses_dir = config['ramses_output_dir']
        output_num = int(os.path.basename(os.path.normpath(ramses_dir)).split('_')[-1])

        info_file = os.path.join(ramses_dir, f"info_{output_num:05d}.txt")
        has_dust = check_simulation_has_dust(info_file)

        if has_dust:
            n_dust = get_dust_species_count(info_file)
        else:
            n_dust = 1
            print("\nAssuming a single dust species derived from gas density.")

        image_config = config['final_image_rendering']
        views_to_render = image_config['views']
        print(f"\nPreparing to render {len(views_to_render)} views: {list(views_to_render.keys())}")

        # Loop through the dictionary items (view_name, view_details).
        for view_name, view_details in views_to_render.items():
            # Generate the command file for the current view.
            cmd_file, image_output_path =\
            create_imaging_cmd_file(config, output_num, n_dust, view_name, view_details)

            # Clean up the specific view directory from any previous run.
            if os.path.exists(image_output_path):
                shutil.rmtree(image_output_path)
            print(f"  - Cleaned up previous directory: {image_output_path}")

            output_base_path = image_config['output_base_path']
            # log_filepath = os.path.join(output_base_path, f'polaris_render_{view_name}.log')
            # log_filepath = os.path.join(output_base_path, f'polaris_render_{view_name}_{output_num:05d}.log')
            base_filename = f'polaris_render_{view_name}_{output_num:05d}'
            # log_filepath = os.path.join(output_base_path, f'{base_filename}.log')
            log_filepath = os.path.join(output_base_path, f'{base_filename}.inner.log')
            print(f"\nExecuting POLARIS for '{view_name}' view... Log will be saved to: '{log_filepath}'\n")

            command = ['polaris', cmd_file]
            with open(log_filepath, 'w') as log_file:
                process = subprocess.Popen(
                    command,
                    stdout = subprocess.PIPE,
                    stderr = subprocess.STDOUT,
                    text = True,
                    encoding = 'utf-8'
                )
                
                for line in process.stdout:
                   sys.stdout.write(line)
                   sys.stdout.flush()
                   log_file.write(line.rstrip('\r\n') + '\n') 

                process.wait()
                if process.returncode != 0:
                    raise subprocess.CalledProcessError(process.returncode, command)
        
        print(f"All final images rendered successfully!")

    except Exception as e:
        print(f"\nERROR in final imaging pipeline: {e}")

if __name__ == '__main__':
    main()

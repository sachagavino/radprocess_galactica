# prepare_radmc3d_inputs.py

# Ilseung Han
#   19.08.2025
#   23.08.2025
#   26.08.2025
#   01.09.2025
#   03.09.2025
#   26.09.2025
#   13.10.2025
#   30.10.2025
#   10.11.2025
#   19.01.2026
#   04.03.2026
#   13.03.2026

import os
import yaml
import numpy as np
import astropy.constants as co

from utils import check_simulation_has_dust, derive_stars_properties, get_dust_species_count

def create_radmc3d_input_files(config: dict, all_stars_derived_props: list, n_dust: int):
    """
    Generates all necessary input files for a RADMC-3D thermal calculation.
    """
    radmc3d_config = config['radmc3d_setup']

    radmc3d_input_dir = radmc3d_config['radmc3d_input_dir']
    print(f"\nGenerating all RADMC-3D input files in: '{radmc3d_input_dir}'")
    os.makedirs(radmc3d_input_dir, exist_ok = True)

    # Generate dustkappa_*.inp.
    # print("  - Converting POLARIS dust opacities...")
    print("  - Creating: dustkappa_*.inp")
    polaris_data_dir = os.path.join(radmc3d_config['polaris_data_input_dir'], 'data')

    kappa_filenames = []
    for i in range(n_dust):
        model_name = f'dust_mixture_{i + 1:03d}'
        kappa_filename = f'dustkappa_{model_name}.inp'
        kappa_filenames.append(model_name)

        polaris_opacity_file = os.path.join(polaris_data_dir, f'{model_name}.dat')
        try:
            # data = np.loadtxt(polaris_opacity_file, skiprows = 27)
            data = np.loadtxt(polaris_opacity_file, skiprows = 29)
            with open(os.path.join(radmc3d_input_dir, kappa_filename), 'w') as f:
                f.write('2\n')
                f.write(f'{len(data)}\n')
                for row in data:
                    # f.write(f'{row[0] / 1e-06:e} {row[-4] * 1e+01:e} {row[-2] * 1e+01:e}\n')
                    f.write(f'{row[0] / 1e-06:e} ' +\
                            f'{(row[-4] * 2 + row[-3]) / 3 * 1e+01:e} ' +\
                            f'{(row[-2] * 2 + row[-1]) / 3 * 1e+01:e}\n')
        except FileNotFoundError:
            raise FileNotFoundError(f"Required POLARIS opacity file not found: {polaris_opacity_file}")

    # Generate dustopac.inp.
    print("  - Creating: dustopac.inp")
    with open(os.path.join(radmc3d_input_dir, 'dustopac.inp'), 'w') as f:
        f.write('2\n')
        f.write(f'{len(kappa_filenames)}\n')
        for name in kappa_filenames:
            f.write('-\n')
            f.write('1\n')
            f.write('0\n')
            f.write(f'{name}\n')

    # Generate wavelength_micron.inp.
    print("  - Creating: wavelength_micron.inp")
    wave_min = radmc3d_config.get('wavelength_min_micron', 0.27)
    wave_max = radmc3d_config.get('wavelength_max_micron', 3000)
    n_wave = radmc3d_config.get('n_wavelengths', 200)
    wavelengths = np.logspace(np.log10(wave_min), np.log10(wave_max), n_wave)

    with open(os.path.join(radmc3d_input_dir, 'wavelength_micron.inp'), 'w') as f:
        f.write(f'{len(wavelengths)}\n')
        for w in wavelengths:
            f.write(f'{w:e}\n')

    # Generate stars.inp
    print("  - Creating: stars.inp")
    with open(os.path.join(radmc3d_input_dir, 'stars.inp'), 'w') as f:
        f.write('2\n')
        f.write(f'{len(all_stars_derived_props)} {len(wavelengths)}\n')

        for star in all_stars_derived_props:
            r_star_cm = star['radius_rsun'] * co.R_sun.cgs.value
            m_star_g = star['mass_msun'] * co.M_sun.cgs.value
            pos_cm = [p * 1e+02 for p in star['pos_m']]
            f.write(f'{r_star_cm: .6e} {m_star_g: .6e} {pos_cm[0]: .6e} {pos_cm[1]: .6e} {pos_cm[2]: .6e}\n')
        
        for w in wavelengths:
            f.write(f'{w:e}\n')

        for star in all_stars_derived_props:
            f.write(f'-{star["temperature_k"]:e}\n')

    # Generate radmc3d.inp.
    print("  - Creating: radmc3d.inp")
    with open(os.path.join(radmc3d_input_dir, 'radmc3d.inp'), 'w') as f:
        n_photons = radmc3d_config.get('n_photons_mctherm', 1000000)
        n_threads = radmc3d_config.get('setthreads', 4)
        f.write(f'nphot = {int(n_photons)}\n')
        f.write(f'nphot_scat = {int(n_photons)}\n')
        f.write(f'setthreads = {n_threads}\n')
        f.write('scattering_mode = 1\n')
        f.write('scattering_mode_max = 1\n')
        # f.write('scattering_mode = 5\n')
        # f.write('scattering_mode_max = 5\n')
        f.write('modified_random_walk = 1\n')
        f.write('rto_style = 3\n')
        f.write('rto_single = 1\n')

def main():
    """Main function to prepare all RADMC-3D input files."""
    config_file_path = 'config.yaml'
    try:
        with open(config_file_path, 'r') as f:
            config = yaml.safe_load(f)
        print(f"Loaded configuration from: '{config_file_path}'")

        ramses_dir = config['ramses_output_dir']
        output_num = int(os.path.basename(os.path.normpath(ramses_dir)).split('_')[-1])
        print(f"Processing RAMSES output number: {output_num}")

        # all_stars = derive_stars_properties(ramses_dir)
        all_stars = derive_stars_properties('/data/pebbles/scratch/ilseung/pipeline/outputs/12-output_00013/')
        if not all_stars: raise ValueError("Could not get star position.")

        info_file = os.path.join(ramses_dir, f"info_{output_num:05d}.txt")
        has_dust = check_simulation_has_dust(info_file)

        if has_dust:
            n_dust = get_dust_species_count(info_file)
        else:
            n_dust = 1
            print("\nAssuming a single dust species derived from gas density.")

        create_radmc3d_input_files(config, all_stars, n_dust)

        print("\nAll RADMC-3D input files have been prepared successfully!")

    except Exception as e:
        print(f"\nERROR in pipeline: {e}")
        raise

if __name__ == '__main__':
    main()

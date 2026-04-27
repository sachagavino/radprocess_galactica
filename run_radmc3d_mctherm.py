# run_radmc3d_mctherm.py

# Ilseung Han
#   21.08.2025
#   23.08.2025
#   13.01.2026
#   19.01.2026

import os
import sys
import yaml
import subprocess

def main():
    """
    Main function to run the RADMC-3D temperature calculation (mctherm).
    """
    config_file_path = 'config.yaml'
    original_directory = os.getcwd() # Save the original directory path.
    log_filepath = '' # Initialize log file path for the error message.

    try:
        with open(config_file_path, 'r') as f:
            config = yaml.safe_load(f)
            print(f"\nLoaded configuration from: '{config_file_path}'")

        radmc3d_run_config = config['radmc3d_run']
        run_directory = radmc3d_run_config['run_directory']

        ramses_dir = config['ramses_output_dir']
        output_num = int(os.path.basename(os.path.normpath(ramses_dir)).split('_')[-1])

        print(f"\nChanging working directory to: '{run_directory}'")
        if not os.path.isdir(run_directory):
            raise FileNotFoundError(f"RADMC-3D run directory not found: '{run_directory}'")
        os.chdir(run_directory)

        # Run RADMC-3D and log the output.
        log_filename = radmc3d_run_config.get('log_filename', f'radmc3d_mctherm_{output_num:05d}.log')
        log_filepath = log_filename

        print(f"\nExecuting 'radmc3d mctherm'... Log will be saved to: '{log_filepath}'")

        command = ['radmc3d', 'mctherm']
        # command = ['/data/pebbles/software/radmc3d-2.0/src/radmc3d', 'mctherm']
        process = subprocess.Popen(
            command,
            stdout = subprocess.PIPE,
            stderr = subprocess.STDOUT,
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

        print("\nRADMC-3D mctherm run completed successfully!")

    except (ValueError, FileNotFoundError, KeyError) as e:
        print(f"\nERROR: Halting pipeline due to a configuration or file issue.")
        print(f"   Details: {e}")
    except subprocess.CalledProcessError as e:
        print(f"\nERROR: RADMC-3D exited with a non-zero status code: {e.returncode}")
        print(f"   Check the log file for details: {os.path.join(run_directory, log_filepath)}")
    except Exception as e:
        print(f"\nERROR: An unexpected error occurred in the pipeline: {e}")
    
    finally:
        # This block will always run, even if an error occurs.
        # It ensures we return to the original directory.
        print(f"\nReturning to original directory: {original_directory}")
        os.chdir(original_directory)

if __name__ == '__main__':
    main()

# update_pymsesrc.py

# Ilseung Han
#   17.08.2025
#   24.08.2025
#   31.10.2025
#   19.01.2026

import os
import json
import yaml

def update_pymsesrc(hydro_file_path: str, pymses_file_path: str):
    """
    Updates the pymsesrc file based on the hydro descriptor.
    If dust_ratio fields are found, it adds them.
    If no dust_ratio fields are found, it removes any existing ones.

    Args:
        hydro_file_path (str): The path to the hydro_file_descriptor.txt file.
        pymses_file_path (str): The path to the pymsesrc file to be updated.
    """
    print(f"\nStarting the '{pymses_file_path}' update process...")
    print(f"\nReading '{os.path.basename(hydro_file_path)}' from:" +\
          f"\n  '{os.path.dirname(hydro_file_path)}'")

    new_dust_fields = []
    try:
        # Open and read the hydro descriptor file line by line.
        with open(hydro_file_path, 'r') as f:
            for line in f:
                # Check if a line contrains the 'dust_ratio' string.
                if 'dust_ratio' in line:
                    parts = line.split(':')
                    var_info, var_name = parts[0], parts[1].strip()
                    # Extract the variable number and convert it to a zero-based index for pymses.
                    ivar = int(var_info.split('#')[1]) - 1
                    # Create the dictionary for the new field and add it to our list.
                    new_dust_fields.append({
                        '__type__': 'scalar_field',
                        '__file_type__': 'hydro',
                        'name': var_name,
                        'ivar': ivar
                    })
                    print(f"  - Found: {var_name} (ivar = {ivar})")

    except FileNotFoundError:
        print(f"\nERROR: Cannot find hydro file at '{hydro_file_path}'")

    try:
        print(f"\nUpdating pymsesrc: '{pymses_file_path}'")
        with open(pymses_file_path, 'r') as f:
            config = json.load(f)

        # Get the original list of AMR fields.
        original_fields = config['RAMSES']['amr_field_descr']

        # Create a new list, excluding any existing dust_ratio fields to avoid duplicates.
        fields_without_dust = [field for field in original_fields
                               if not field.get('name', '').startswith('dust_ratio')]

        # Rebuild the final list of fields.
        final_amr_fields = []
        for field in fields_without_dust:
            # Add the current field from the original list.
            final_amr_fields.append(field)
            # If the field we just added was the pressure ('P') field...
            if field.get('name') == 'P':
                # ...insert all the new dust fields we found after it.
                final_amr_fields.extend(new_dust_fields)
                # print("- Inserting dust fields after 'P' field...")

        # Replace the old field list in the configuration with the our newly constructed list.
        config['RAMSES']['amr_field_descr'] = final_amr_fields

        # Write the updated configuration dictionary back to the pymsesrc file.
        with open(pymses_file_path, 'w') as f:
            json.dump(config, f, indent = 4)
        print(f"\nSuccessfully updated '{pymses_file_path}'")

    except Exception as e:
        print(f"\nERROR: An unexpected error occurred while updated '{pymses_file_path}': {e}")
        raise

def main():
    """
    Main function to update the pymsesrc file.
    """
    config_file_path = 'config.yaml'
    try:
        with open(config_file_path, 'r') as f:
            config = yaml.safe_load(f)
            print(f"\nLoaded configuration from: '{config_file_path}'")

        ramses_dir = config['ramses_output_dir']
        hydro_file = os.path.join(ramses_dir, 'hydro_file_descriptor.txt')
        pymses_file = config['pymsesrc_path']

        update_pymsesrc(hydro_file, pymses_file)

    except FileNotFoundError:
        print(f"\nERROR: Configuration file not found at '{config_file_path}'")
        raise
    except Exception as e:
        print(f"\nERROR: Failed to process configuration: {e}")
        raise

if __name__ == '__main__':
    main()

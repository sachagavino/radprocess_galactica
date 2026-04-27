import subprocess

for script in [
    "update_pymsesrc.py",
    "create_polaris_grid.py",
    "create_radmc3d_grid.py",
    "run_polaris_opacity.py",
    "prepare_radmc3d_inputs.py",
    "run_radmc3d_mctherm.py",
    "merge_temperature.py",
    "render_final_images.py",
]:
   subprocess.run(["python", script], check=True)
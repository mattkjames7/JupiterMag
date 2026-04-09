import numpy as np
import JupiterMag as jm
import json
import sys
import os

sys.path.append("..")
from common import get_con2020_field


# 10 cartesian points in System III
x = np.array([0, 5.0, 10.0, 20.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0, 50.0, 100.0, 50.0, 100.0])
y = np.array([5.0, 0, 0, 0, 5.0, 10.0, 20.0, 10.0, 5.0, 10.0, 50.0, 100.0, 100.0, 50.0])
z = np.array([20.0, 0, 0, 0, 0, 0, 0, 5.0, 10.0, 20.0, 1.0, 1.0, 20.0, 20.0])

# Spherical points in System III
r = np.sqrt(x**2 + y**2 + z**2)
theta = np.arccos(z/r)
phi = np.arctan2(y,x)


cfg_vars = {
    'mu_i': [139.6, 130.0, 150.0],
    'i_rho': [16.7, 0.0, 30.0],
    'r0': [7.8, 5.0, 10.0],
    'r1': [51.4, 30.0, 70.0],
    'd': [3.6, 2.0, 5.0],
    'xt': [9.3, 0.0, 20.0],
    'xp': [155.8, 90.0, 180.0],
    'Edwards': [True, False],
    'error_check': [True],
    'CartesianIn': [True, False],
    'CartesianOut': [True, False],
    'equation_type': ['hybrid', "analytic", "integral"],
    'Smooth': [False, True],
    'DeltaRho': [1.0],
    'DeltaZ': [0.1],
    'g': [417659.3836476442],
    'azfunc': ['connerney', "lmic"],
    'wO_open': [0.1],
    'wO_om': [0.35],
    'thetamm': [16.1],
    'dthetamm': [0.5],
    'thetaoc': [10.716],
    'dthetaoc': [0.125]
}


def generate_config_list():

    out = []

    # list all combinations of Edwards, CartesianIn, CartesianOut, azfunc, and Smooth
    full_keys = ['Edwards', 'CartesianIn', 'CartesianOut', 'azfunc', 'Smooth']
    full_combinations = []
    for Edwards in cfg_vars['Edwards']:
        for CartesianIn in cfg_vars['CartesianIn']:
            for CartesianOut in cfg_vars['CartesianOut']:
                for azfunc in cfg_vars['azfunc']:
                    for Smooth in cfg_vars['Smooth']:
                        full_combinations.append({
                            'Edwards': Edwards,
                            'CartesianIn': CartesianIn,
                            'CartesianOut': CartesianOut,
                            'azfunc': azfunc,
                            'Smooth': Smooth
                        })

    # these keys will be tested for each of their values, but with the first element of the opther part keys (but also for every combination of the full keys)
    part_keys = [key for key in cfg_vars.keys() if key not in full_keys]

    for combination in full_combinations:
        for key in part_keys:
            for value in cfg_vars[key]:
                config = combination.copy()
                config[key] = value
                for k, v in cfg_vars.items():
                    if k not in config:
                        config[k] = v[0]
                out.append(config)
    return out


def save_con2020_data(filename, overwrite=False):
    if not overwrite and os.path.exists(filename):
        print(f"File {filename} already exists. Use --overwrite to overwrite it.")
        return

    print(f"Saving Con2020 field test data to {filename}...")

    configs = generate_config_list()

    data = []
    for i, cfg in enumerate(configs):
        print(f"\rProcessing config {i+1}/{len(configs)}", end="")
        field = get_con2020_field(x, y, z, cfg)
        data.append({
            "input": {
                "args": [x.tolist(), y.tolist(), z.tolist(), cfg],
            },
            "output": {
                "result": [field[0].tolist(), field[1].tolist(), field[2].tolist()]
            }
        })
    print()  # Move to the next line after the progress bar

    with open(filename, "w") as f:
        json.dump(data, f, indent=2)

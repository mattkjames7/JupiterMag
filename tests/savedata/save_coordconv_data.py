import numpy as np
import JupiterMag as jm
import json
import os

# 10 cartesian points in System III
x = np.array([0, 5.0, 10.0, 20.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0])
y = np.array([5.0, 0, 0, 0, 5.0, 10.0, 20.0, 10.0, 5.0, 10.0])
z = np.array([20.0, 0, 0, 0, 0, 0, 0, 5.0, 10.0, 20.0])


# some different values for xt and xp
xt = np.array([0, 5.0, 10.0, 20.0, 45.0, 90.0, 135.0])
xp = np.array([0, 5.0, 10.0, 20.0, 45.0, 90.0, 135.0, 180.0, 270.0])


def save_coordconv_data(filename, overwrite=False):
    if not overwrite and os.path.exists(filename):
        print(f"File {filename} already exists. Use --overwrite to overwrite it.")
        return

    print(f"Saving coordinate conversion test data to {filename}...")

    functions = {
        "MagtoSIII": jm.CoordConv.MagtoSIII,
        "SIIItoMag": jm.CoordConv.SIIItoMag,
    }

    data = []
    for func_name, func in functions.items():
        for xt_val in xt:
            for xp_val in xp:
                ox, oy, oz = func(x, y, z, xt_val, xp_val)
                for i in range(len(x)):
                    data.append(
                        {
                            "function": func_name,
                            "input": {
                                "args": [x[i], y[i], z[i], xt_val, xp_val],
                                "kwargs": {},
                            },
                            "output": {"result": [ox[i], oy[i], oz[i]]},
                        }
                    )

    with open(filename, "w") as f:
        json.dump(data, f, indent=2)

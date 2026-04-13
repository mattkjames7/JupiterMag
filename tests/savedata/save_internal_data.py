import numpy as np
import JupiterMag as jm
import json
import sys
import os

sys.path.append("..")
from common import get_internal_field


# 10 cartesian points in System III
x = np.array([0, 5.0, 10.0, 20.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0])
y = np.array([5.0, 0, 0, 0, 5.0, 10.0, 20.0, 10.0, 5.0, 10.0])
z = np.array([20.0, 0, 0, 0, 0, 0, 0, 5.0, 10.0, 20.0])

# Spherical points in System III
r = np.sqrt(x**2 + y**2 + z**2)
theta = np.arccos(z/r)
phi = np.arctan2(y,x)


model_configs = {
    "jrm09": [
        {
            "CartIn": True,
            "CartOut": True,
            "Degree": 10
        },
        {
            "CartIn": True,
            "CartOut": True,
            "Degree": 5
        },
        {
            "CartIn": True,
            "CartOut": True,
            "Degree": 1
        },
        {
            "CartIn": True,
            "CartOut": False,
            "Degree": 10
        },
        {
            "CartIn": True,
            "CartOut": False,
            "Degree": 5
        },
        {
            "CartIn": True,
            "CartOut": False,
            "Degree": 1
        },
        {
            "CartIn": False,
            "CartOut": True,
            "Degree": 10
        },
        {
            "CartIn": False,
            "CartOut": True,
            "Degree": 5
        },
        {
            "CartIn": False,
            "CartOut": True,
            "Degree": 1
        },
        {
            "CartIn": False,
            "CartOut": False,
            "Degree": 10
        },
        {
            "CartIn": False,
            "CartOut": False,
            "Degree": 5
        },
        {
            "CartIn": False,
            "CartOut": False,
            "Degree": 1
        }
    ],
    "vip4": [
        {
            "CartIn": True,
            "CartOut": True,
            "Degree": 4
        },
        {
            "CartIn": True,
            "CartOut": False,
            "Degree": 4
        },
        {
            "CartIn": False,
            "CartOut": True,
            "Degree": 4
        },
        {
            "CartIn": False,
            "CartOut": False,
            "Degree": 4
        }
    ]
}


def save_internal_field_data(filename, overwrite=False):
    if not overwrite and os.path.exists(filename):
        print(f"File {filename} already exists. Use --overwrite to overwrite it.")
        return

    print(f"Saving internal field test data to {filename}...")

    data = []
    for model_name, configs in model_configs.items():
        for cfg in configs:
            field = get_internal_field(x, y, z, model=model_name, CartIn=cfg["CartIn"], CartOut=cfg["CartOut"], MaxDeg=cfg["Degree"])
            data.append({
                "function": "Internal.Field",
                "input": {
                    "args": [x.tolist(), y.tolist(), z.tolist()],
                    "kwargs": {
                        "model": model_name,
                        "CartIn": cfg["CartIn"],
                        "CartOut": cfg["CartOut"],
                        "MaxDeg": cfg["Degree"]
                    }
                },
                "output": {
                    "result": [item.tolist() for item in field]
                }
            })

    with open(filename, "w") as f:
        json.dump(data, f, indent=2)

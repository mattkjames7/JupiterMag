import numpy as np
import JupiterMag as jm
import json
import sys
import os

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(CURRENT_DIR, ".."))
from common import get_trace_footprints


internal_models = ["jrm09", "vip4"]
external_models = ["Con2020", "none"]  # TODO: make model names case-insensitive
equation_types = ["analytic", "hybrid", "integral"]  # integral is too slow for testing

# start positions in System III
r = np.array([2.0, 5.0, 10.0, 20.0, 50.0, 60.0])
z = np.array([0.0, 0.0, 5.0, 10.0, 5.0, -10.0])
phi = np.array([0.0, 45.0, 90.0, 135.0, 180.0, 270.0])

x = r * np.cos(np.radians(phi))
y = r * np.sin(np.radians(phi))


def generate_inputs():
    configs = []
    for int_model in internal_models:
        for ext_model in external_models:
            if ext_model == "none":
                eq_types = ["analytic"]  # only test analytic for no external field
            else:
                eq_types = equation_types
            for eq_type in eq_types:
                configs.append({
                    "IntModel": int_model,
                    "ExtModel": ext_model,
                    "equation_type": eq_type
                })

    test_input_data = []
    for cfg in configs:
        for i in range(len(x)):

            # skip the integral cases for r > 2 to save time, since they are very slow
            if cfg["equation_type"] == "integral" and r[i] > 2.0:
                continue

            test_input_data.append({
                "function": "TraceFootprints",
                "input": {
                    "args": [x[i], y[i], z[i], cfg["IntModel"], cfg["ExtModel"], {"equation_type": cfg["equation_type"]}],
                    "kwargs": {}
                },
                "output": {
                    "result": None  # to be filled in with expected result
                }
            })

    return test_input_data


def save_trace_data(filename, overwrite=False):
    if not overwrite and os.path.exists(filename):
        print(f"File {filename} already exists. Use --overwrite to overwrite it.")
        return

    print(f"Saving TraceFootprints test data to {filename}...")

    test_input_data = generate_inputs()

    for test_case in test_input_data:

        # calculate r
        r = float(np.sqrt(test_case["input"]["args"][0]**2 + test_case["input"]["args"][1]**2 + test_case["input"]["args"][2]**2))
        r_str = f"{r:5.1f}"

        # get equation_type
        eq_type = test_case["input"]["args"][5]["equation_type"]

        print(f"\rProcessing config {test_input_data.index(test_case)+1}/{len(test_input_data)} r={r_str} equation_type={eq_type:<8}", end="")
        args = test_case["input"]["args"]
        kwargs = test_case["input"]["kwargs"]
        result = get_trace_footprints(*args, **kwargs)
        test_case["output"]["result"] = result
    print("\nDone.")

    with open(filename, "w") as f:
        json.dump(test_input_data, f, indent=2)

import JupiterMag as jm
import os
import json
import numpy as np
import argparse
from save_modelfield_data import save_modelfield_data
from save_coordconv_data import save_coordconv_data
from save_internal_data import save_internal_field_data
from save_con2020_data import save_con2020_data

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_DIR = os.path.join(SCRIPT_DIR, "../data")


"""
Test data format:
{
    "function": "function_name",
    "input": {
        "args": [],
        "kwargs": {}
    },
    "output": {
        "result": "expected result"
    }
}

"""


def main():
    parser = argparse.ArgumentParser(description="Generate test data for JupiterMag functions")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing test data files")
    args = parser.parse_args()


    save_modelfield_data(os.path.join(TEST_DATA_DIR, "modelfield_data.json"), overwrite=args.overwrite)
    save_coordconv_data(os.path.join(TEST_DATA_DIR, "coordconv_data.json"), overwrite=args.overwrite)
    save_internal_field_data(os.path.join(TEST_DATA_DIR, "internal_field_data.json"), overwrite=args.overwrite)
    save_con2020_data(os.path.join(TEST_DATA_DIR, "con2020_field_data.json"), overwrite=args.overwrite)


if __name__ == "__main__":
    main()

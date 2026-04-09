import JupiterMag as jm
import os
import json
import numpy as np
from save_modelfield_data import save_modelfield_data
from save_coordconv_data import save_coordconv_data
from save_internal_data import save_internal_field_data

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
    save_modelfield_data(os.path.join(TEST_DATA_DIR, "modelfield_data.json"))
    save_coordconv_data(os.path.join(TEST_DATA_DIR, "coordconv_data.json"))
    save_internal_field_data(os.path.join(TEST_DATA_DIR, "internal_field_data.json"))


if __name__ == "__main__":
    main()

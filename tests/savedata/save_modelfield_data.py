import numpy as np
import json
import sys
import os

sys.path.append("..")
from common import get_model_field


# 10 cartesian points in System III
x = np.array([0, 5.0, 10.0, 20.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0])
y = np.array([5.0, 0, 0, 0, 5.0, 10.0, 20.0, 10.0, 5.0, 10.0])
z = np.array([20.0, 0, 0, 0, 0, 0, 0, 5.0, 10.0, 20.0])

# Spherical points in System III
r = np.sqrt(x**2 + y**2 + z**2)
theta = np.arccos(z/r)
phi = np.arctan2(y,x)

test_input_data = [
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": "Con2020", "CartIn": True, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": "Con2020", "CartIn": True, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": "Con2020", "CartIn": True, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": "Con2020", "CartIn": True, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": None, "CartIn": True, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": None, "CartIn": True, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": None, "CartIn": True, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": None, "CartIn": True, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": "Con2020", "CartIn": False, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": "Con2020", "CartIn": False, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": "Con2020", "CartIn": False, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": "Con2020", "CartIn": False, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": "none", "CartIn": False, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": "none", "CartIn": False, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": "none", "CartIn": False, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": "none", "CartIn": False, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    }
]


def save_modelfield_data(filename, overwrite=False):
    if not overwrite and os.path.exists(filename):
        print(f"File {filename} already exists. Use --overwrite to overwrite it.")
        return

    print(f"Saving ModelField test data to {filename}...")

    for test_case in test_input_data:
        args = test_case["input"]["args"]
        kwargs = test_case["input"]["kwargs"]

        bx, by, bz = get_model_field(*args, **kwargs)
        test_case["output"]["result"] = [bx.tolist(), by.tolist(), bz.tolist()]  # convert to list for JSON serialization

    with open(filename, "w") as f:
        json.dump(test_input_data, f, indent=2)

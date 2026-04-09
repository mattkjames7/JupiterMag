import numpy as np
import JupiterMag as jm
import json


# 10 cartesian points in System III
x = np.array([0, 5.0, 10.0, 20.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0])
y = np.array([5.0, 0, 0, 0, 5.0, 10.0, 20.0, 10.0, 5.0, 10.0])
z = np.array([20.0, 0, 0, 0, 0, 0, 0, 5.0, 10.0, 20.0])

# Spherical points in System III
r = np.sqrt(x**2 + y**2 + z**2)
theta = np.arccos(z/r)
phi = np.arctan2(y,x)

# Default configs for the models we will be testing
default_jrm09_cfg = jm.Internal.Config("default", Model="jrm09")
default_vip4_cfg = jm.Internal.Config("default", Model="vip4")
default_con2020_cfg = jm.Con2020.Config("default")

test_input_data = [
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": "con2020", "CartIn": True, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": "con2020", "CartIn": True, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": "con2020", "CartIn": True, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [x.tolist(), y.tolist(), z.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": "con2020", "CartIn": True, "CartOut": False}
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
            "kwargs": {"IntModel": "jrm09", "ExtModel": "con2020", "CartIn": False, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "jrm09", "ExtModel": "con2020", "CartIn": False, "CartOut": True}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": "con2020", "CartIn": False, "CartOut": False}
        },
        "output": {
            "result": None  # to be filled in with expected result
        }
    },
    {
        "function": "ModelField",
        "input": {
            "args": [r.tolist(), theta.tolist(), phi.tolist()],
            "kwargs": {"IntModel": "vip4", "ExtModel": "con2020", "CartIn": False, "CartOut": True}
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


# safe model field function wrapper
def get_model_field(x, y, z, IntModel="jrm09", ExtModel="con2020", CartIn=True, CartOut=True):
    
    if IntModel == "jrm09":
        IntModel_cfg = default_jrm09_cfg
    elif IntModel == "vip4":
        IntModel_cfg = default_vip4_cfg
    else:
        raise ValueError(f"Invalid internal model: {IntModel}")
    
    jm.Internal.Config(**IntModel_cfg)
    jm.Con2020.Config(**default_con2020_cfg)

    return jm.ModelField(x, y, z, IntModel=IntModel, ExtModel=ExtModel, CartIn=CartIn, CartOut=CartOut)


def save_modelfield_data(filename):
    for test_case in test_input_data:
        args = test_case["input"]["args"]
        kwargs = test_case["input"]["kwargs"]

        bx, by, bz = get_model_field(*args, **kwargs)
        test_case["output"]["result"] = [bx.tolist(), by.tolist(), bz.tolist()]  # convert to list for JSON serialization

    with open(filename, "w") as f:
        json.dump(test_input_data, f, indent=2)
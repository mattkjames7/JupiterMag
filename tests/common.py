import JupiterMag as jm


# Default configs for the models we will be testing
default_jrm09_cfg = jm.Internal.Config("default", Model="jrm09")
default_vip4_cfg = jm.Internal.Config("default", Model="vip4")
default_con2020_cfg = jm.Con2020.Config("default")


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


def get_internal_field(p0, p1, p2, model, CartIn=True, CartOut=True, MaxDeg=None):

    cfg = {
        "Model": model,
        "CartesianIn": CartIn,
        "CartesianOut": CartOut,
        "Degree": MaxDeg if MaxDeg is not None else 0  # 0 will be treated as default by the C function
    }
    jm.Internal.Config(**cfg)

    return jm.Internal.Field(p0, p1, p2)


def get_con2020_field(p0, p1, p2, cfg):

    jm.Con2020.Config(**cfg)

    return jm.Con2020.Field(p0, p1, p2)

# TODO: this looks like something isn't actually working correctly
def get_trace_footprints(x, y, z, internal_model, external_model, con2020_cfg):

    jm.Internal.Config("default")
    jm.Internal.Config(Model=internal_model)

    jm.Con2020.Config("default")
    jm.Con2020.Config(**con2020_cfg)

    T = jm.TraceField(x, y, z, IntModel=internal_model, ExtModel=external_model)

    fp_arrays = {
        "ionosphere": T.ionosphere,
        "surface": T.surface,
        "equator": T.equator
    }

    out = {}
    for key, arr in fp_arrays.items():
        out[key] = {}
        names = arr.dtype.names
        for i, name in enumerate(names):
            out[key][name] = arr[name].tolist()

    return out

import JupiterMag as jm

# Default configs are computed lazily to avoid mutating global C++
# configuration state when this module is imported.
_default_internal_cfgs = {}
_default_con2020_cfg = None


def _get_default_internal_cfg(model):
    cfg = _default_internal_cfgs.get(model)
    if cfg is None:
        current_cfg = jm.Internal.Config()
        try:
            cfg = dict(jm.Internal.Config("default", Model=model))
        finally:
            jm.Internal.Config(**current_cfg)
        _default_internal_cfgs[model] = cfg

    return dict(_default_internal_cfgs[model])


def _get_default_con2020_cfg():
    global _default_con2020_cfg

    if _default_con2020_cfg is None:
        current_cfg = jm.Con2020.Config()
        try:
            _default_con2020_cfg = dict(jm.Con2020.Config("default"))
        finally:
            jm.Con2020.Config(**current_cfg)

    return dict(_default_con2020_cfg)


# safe model field function wrapper
def get_model_field(x, y, z, IntModel="jrm09", ExtModel="con2020", CartIn=True, CartOut=True):

    if IntModel == "jrm09":
        IntModel_cfg = _get_default_internal_cfg("jrm09")
    elif IntModel == "vip4":
        IntModel_cfg = _get_default_internal_cfg("vip4")
    else:
        raise ValueError(f"Invalid internal model: {IntModel}")

    current_internal_cfg = jm.Internal.Config()
    current_con2020_cfg = jm.Con2020.Config()
    try:
        jm.Internal.Config(**IntModel_cfg)
        jm.Con2020.Config(**_get_default_con2020_cfg())

        return jm.ModelField(
            x,
            y,
            z,
            IntModel=IntModel,
            ExtModel=ExtModel,
            CartIn=CartIn,
            CartOut=CartOut,
        )
    finally:
        jm.Internal.Config(**current_internal_cfg)
        jm.Con2020.Config(**current_con2020_cfg)


def get_internal_field(p0, p1, p2, model, CartIn=True, CartOut=True, MaxDeg=None):

    cfg = {
        "Model": model,
        "CartesianIn": CartIn,
        "CartesianOut": CartOut,
        "Degree": (MaxDeg if MaxDeg is not None else 0),  # 0 will be treated as default by the C function
    }
    current_internal_cfg = jm.Internal.Config()
    try:
        jm.Internal.Config(**cfg)
        return jm.Internal.Field(p0, p1, p2)
    finally:
        jm.Internal.Config(**current_internal_cfg)


def get_con2020_field(p0, p1, p2, cfg):

    current_con2020_cfg = jm.Con2020.Config()
    try:
        jm.Con2020.Config(**cfg)
        return jm.Con2020.Field(p0, p1, p2)
    finally:
        jm.Con2020.Config(**current_con2020_cfg)


def get_trace_footprints(x, y, z, internal_model, external_model, con2020_cfg):

    current_internal_cfg = jm.Internal.Config()
    current_con2020_cfg = jm.Con2020.Config()
    try:
        jm.Internal.Config("default")
        jm.Internal.Config(Model=internal_model)

        jm.Con2020.Config("default")
        jm.Con2020.Config(**con2020_cfg)

        T = jm.TraceField(x, y, z, IntModel=internal_model, ExtModel=external_model)

        fp_arrays = {
            "ionosphere": T.ionosphere,
            "surface": T.surface,
            "equator": T.equator,
        }

        out = {}
        for key, arr in fp_arrays.items():
            out[key] = {}
            names = arr.dtype.names
            for i, name in enumerate(names):
                out[key][name] = arr[name].tolist()

        return out
    finally:
        jm.Internal.Config(**current_internal_cfg)
        jm.Con2020.Config(**current_con2020_cfg)

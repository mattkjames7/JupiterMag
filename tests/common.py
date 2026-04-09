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

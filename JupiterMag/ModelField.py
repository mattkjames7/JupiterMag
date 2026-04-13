import numpy as np
from ._CFunctions import _CModelFieldArray


def ModelField(x, y, z, IntModel="jrm09", ExtModel="con2020", CartIn=True, CartOut=True):
    '''
    Calculate the magnetic field at a given location using the internal and
    external models specified. The internal and external models are calculated
    separately and then added together to get the total field.
    
    Inputs
    ======
    x : float or array-like
    x coordinate(s) in System III (km)
    y : float or array-like
    y coordinate(s) in System III (km)
    z : float or array-like
    z coordinate(s) in System III (km)
    IntModel : str, optional
    Name of internal model to use. Default is "jrm09".
    ExtModel : str, optional
    Name of external model to use. Default is "con2020".
    CartIn : bool, optional
    If True, input coordinates are in Cartesian System III. If False, input 
    coordinates are in spherical System III (r, theta, phi). Default is True.
    CartOut : bool, optional
    If True, output field components are in Cartesian System III. If False, output 
    field components are in spherical System III (Br, Bt, Bp). Default is True.
    
    Returns
    =======
    Bx : float or array-like
    x component of magnetic field in System III (nT)
    By : float or array-like
    y component of magnetic field in System III (nT)
    Bz : float or array-like
    z component of magnetic field in System III (nT)
    
    '''
    
    x = np.atleast_1d(x).astype(np.float64)
    y = np.atleast_1d(y).astype(np.float64)
    z = np.atleast_1d(z).astype(np.float64)

    # Ensure all inputs are shape-compatible and contiguous before passing
    # pointers and length to the C layer.
    x, y, z = np.broadcast_arrays(x, y, z)
    x = np.ascontiguousarray(x, dtype=np.float64).ravel()
    y = np.ascontiguousarray(y, dtype=np.float64).ravel()
    z = np.ascontiguousarray(z, dtype=np.float64).ravel()
    n = np.int32(x.size)
    
    Bx = np.zeros(n,dtype=np.float64)
    By = np.zeros(n,dtype=np.float64)
    Bz = np.zeros(n,dtype=np.float64)

    IntModel_c = "none".encode('utf-8') if IntModel is None else IntModel.encode('utf-8')

    if ExtModel is None:
        ExtModel_c = "none".encode('utf-8')
    elif ExtModel.lower() == "con2020":
        ExtModel_c = "Con2020".encode('utf-8')
    else:
        ExtModel_c = ExtModel.encode('utf-8')
    
    _CModelFieldArray(
        n,x,y,z,
        IntModel_c,
        ExtModel_c,
        CartIn,
        CartOut,
        Bx,
        By,
        Bz
    )
    
    return Bx, By, Bz

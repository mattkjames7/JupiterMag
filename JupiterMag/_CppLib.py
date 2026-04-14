import numpy as np
import os
import subprocess
import ctypes
import platform
import fnmatch
from . import Globals


def _LibPath():
    """
    Return a path to the C++ library

    Returns
    =======
    path : str
            path to the library file.

    """
    return os.path.dirname(__file__) + "/__data/libjupitermag/lib/"


def _LibPaths():
    """
    Return candidate directories containing the native library.

    Returns
    =======
    paths : list
            Candidate directories.

    """
    root = os.path.dirname(__file__) + "/__data/libjupitermag/"
    return [root + "lib/", root + "lib64/", root + "bin/"]


def _LibName(WithPath=False):
    """
    Return the name of the library.

    Inputs
    ======
    WithPath : bool
            If True then the full path to the library will be included.

    Returns
    =======
    libpath : str
            Library name.

    """
    path = _LibPath() if WithPath else ""

    osname = platform.uname().system
    libexts = {"Linux": "so", "Windows": "dll", "Darwin": "dylib"}

    ext = libexts[osname]

    if ext is None:
        raise Exception("The Operating System ({:s}) is not supported".format(osname))

    return path + "libjupitermag." + ext


def _LibNameCandidates(WithPath=False):
    """
    Return one or more candidate library names.

    Inputs
    ======
    WithPath : bool
            If True then full paths to candidate libraries are returned.

    Returns
    =======
    out : list
            Candidate library names or full paths.

    """
    osname = platform.uname().system
    libexts = {"Linux": "so", "Windows": "dll", "Darwin": "dylib"}

    ext = libexts[osname]

    if ext is None:
        raise Exception("The Operating System ({:s}) is not supported".format(osname))

    name = "libjupitermag." + ext
    if not WithPath:
        return [name]

    return [p + name for p in _LibPaths()]


def _LibExists():
    """
    Check if the library file exists.

    Returns
    =======
    exists : bool
            True if the file exists
    """
    for name in _LibNameCandidates(True):
        if os.path.isfile(name):
            return True
    return False


def getWindowsSearchPaths():
    """Scan the directories within PATH and look for std C++ libs"""
    paths = os.getenv("PATH")
    paths = paths.split(";")

    pattern = "libstdc++*.dll"

    out = []
    for p in paths:
        if os.path.isdir(p):
            files = os.listdir(p)
            mch = any(fnmatch.fnmatch(f, pattern) for f in files)
            if mch:
                out.append(p)

    return out


def addWindowsSearchPaths():

    paths = getWindowsSearchPaths()
    for p in paths:
        if os.path.isdir(p):
            os.add_dll_directory(p)


def _GetLib():
    """
    Return an instance of the C++ library

    Returns
    =======
    lib : ctypes.CDLL
            C++ library containing the field model code
    """
    fnames = _LibNameCandidates(True)
    fname = None
    for n in fnames:
        if os.path.isfile(n):
            fname = n
            break
    if fname is None:
        fname = _LibName(True)

    try:
        print("Importing Library")
        if platform.system() == "Windows":
            addWindowsSearchPaths()
            lib = ctypes.CDLL(fname)
        elif platform.system() == "Darwin":
            cwd = os.getcwd()
            loaded = False
            for path in _LibPaths():
                if not os.path.isdir(path):
                    continue
                try:
                    os.chdir(path)
                    lib = ctypes.CDLL(_LibName(False))
                    loaded = True
                    break
                except Exception:
                    continue
            os.chdir(cwd)
            if not loaded:
                raise OSError("Unable to load native JupiterMag library")
        else:
            lib = ctypes.CDLL(fname)
        print("done")
    except Exception as e:
        print("Importing C++ library failed. Please reinstall...")
        print(e)
        raise SystemExit

    return lib

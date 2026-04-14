import os
import ctypes
import platform
import fnmatch

_LIB_EXTENSIONS = {"Linux": "so", "Windows": "dll", "Darwin": "dylib"}


def _lib_extension():
    ext = _LIB_EXTENSIONS.get(platform.system())
    if ext is None:
        raise RuntimeError(f"The Operating System ({platform.system():s}) is not supported")
    return ext


def _lib_filename():
    return f"libjupitermag.{_lib_extension()}"


def _candidate_lib_paths():
    root = os.path.join(os.path.dirname(__file__), "__data", "libjupitermag")
    name = _lib_filename()
    return [
        os.path.join(root, "lib", name),
        os.path.join(root, "lib64", name),
        os.path.join(root, "bin", name),
    ]


def _LibPath():
    """
    Return a path to the C++ library

    Returns
    =======
    path : str
            path to the library file.

    """
    return os.path.join(os.path.dirname(__file__), "__data", "libjupitermag", "lib") + "/"


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
    name = _lib_filename()
    return os.path.join(_LibPath(), name) if WithPath else name


def _LibExists():
    """
    Check if the library file exists.

    Returns
    =======
    exists : bool
            True if the file exists
    """
    return any(os.path.isfile(path) for path in _candidate_lib_paths())


def _add_windows_search_paths():
    """Scan PATH entries and add directories containing libstdc++ DLLs."""
    pattern = "libstdc++*.dll"
    for entry in os.getenv("PATH", "").split(";"):
        if not os.path.isdir(entry):
            continue
        files = os.listdir(entry)
        if any(fnmatch.fnmatch(name, pattern) for name in files):
            os.add_dll_directory(entry)


def _GetLib():
    """
    Return an instance of the C++ library

    Returns
    =======
    lib : ctypes.CDLL
            C++ library containing the field model code
    """
    candidates = _candidate_lib_paths()
    fname = next((path for path in candidates if os.path.isfile(path)), _LibName(True))

    try:
        print("Importing Library")
        if platform.system() == "Windows":
            _add_windows_search_paths()
            lib = ctypes.CDLL(fname)
        elif platform.system() == "Darwin":
            cwd = os.getcwd()
            loaded = False
            for path in candidates:
                lib_dir = os.path.dirname(path)
                if not os.path.isdir(lib_dir):
                    continue
                try:
                    os.chdir(lib_dir)
                    lib = ctypes.CDLL(_LibName(False))
                    loaded = True
                    break
                except Exception:
                    continue
                finally:
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

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from setuptools.dist import Distribution
import subprocess
import os
import platform
import shutil
from pathlib import Path
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel


class bdist_wheel(_bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        # Mark wheel as platform-specific because it contains native binaries.
        self.root_is_pure = False


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        # Force platlib wheel layout for bundled native shared libraries.
        return True


class CustomBuild(build_py):
    def run(self):
        self.execute(self.target_build, ())
        build_py.run(self)
        self.copy_native_artifacts()

    def target_build(self):
        root = Path(__file__).resolve().parent
        lib_root = root / "JupiterMag" / "__data" / "libjupitermag"
        build_cmd = self.get_finalized_command("build")
        build_dir = Path(build_cmd.build_temp) / "libjupitermag-cmake"
        build_dir.mkdir(parents=True, exist_ok=True)
        expected_lib = self._find_main_library(lib_root)

        cmake_configure = [
            "cmake",
            "-S",
            str(lib_root),
            "-B",
            str(build_dir),
            f"-DCMAKE_INSTALL_PREFIX={lib_root}",
            "-DBUILD_SHARED_LIBS=ON",
            "-DLIBJUPITERMAG_BUILD_TESTS=OFF",
        ]

        # Prefer Ninja when available, unless caller has explicitly selected
        # a CMake generator.
        if shutil.which("ninja") and "CMAKE_GENERATOR" not in os.environ:
            cmake_configure += ["-G", "Ninja"]

        subprocess.check_call(cmake_configure, stderr=subprocess.STDOUT)

        cmake_build = ["cmake", "--build", str(build_dir)]
        if platform.system() == "Windows":
            cmake_build += ["--config", "Release"]
        subprocess.check_call(cmake_build, stderr=subprocess.STDOUT)

        cmake_install = ["cmake", "--install", str(build_dir)]
        if platform.system() == "Windows":
            cmake_install += ["--config", "Release"]
        subprocess.check_call(cmake_install, stderr=subprocess.STDOUT)

        if platform.system() == "Windows":
            # CMake installs runtime DLLs under bin/ on Windows; keep a copy in
            # lib/ so existing Python loader logic continues to work unchanged.
            dll_in_bin = lib_root / "bin" / "libjupitermag.dll"
            dll_in_lib = lib_root / "lib" / "libjupitermag.dll"
            if dll_in_bin.is_file() and not dll_in_lib.is_file():
                dll_in_lib.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(dll_in_bin, dll_in_lib)

        # Normalize Linux manylinux installs that may use lib64/ so Python loader
        # and package_data can continue reading from lib/.
        self._ensure_main_library_in_lib(lib_root)
        expected_lib = self._find_main_library(lib_root)

        # Fail fast if the native library was not produced.
        if not expected_lib.is_file():
            raise RuntimeError(f"Native library was not built at {expected_lib}")

    def _main_library_candidates(self, lib_root: Path):
        ext = {
            "Windows": "dll",
            "Linux": "so",
            "Darwin": "dylib",
        }.get(platform.system())
        if ext is None:
            raise RuntimeError(f"Unsupported OS: {platform.system()}")

        name = f"libjupitermag.{ext}"
        return [
            lib_root / "lib" / name,
            lib_root / "lib64" / name,
            lib_root / "bin" / name,
        ]

    def _find_main_library(self, lib_root: Path) -> Path:
        candidates = self._main_library_candidates(lib_root)
        for path in candidates:
            if path.is_file():
                return path
        return candidates[0]

    def _ensure_main_library_in_lib(self, lib_root: Path):
        lib_path = lib_root / "lib"
        preferred = self._main_library_candidates(lib_root)[0]
        if preferred.is_file():
            return

        found = self._find_main_library(lib_root)
        if found.is_file() and found.parent != lib_path:
            lib_path.mkdir(parents=True, exist_ok=True)
            shutil.copy2(found, preferred)

    def copy_native_artifacts(self):
        """Copy generated native artifacts into build_lib for wheel packaging."""
        exts = {".dll", ".so", ".dylib"}
        root = Path(__file__).resolve().parent
        data_root = root / "JupiterMag" / "__data" / "libjupitermag"
        dst_root = Path(self.build_lib) / "JupiterMag" / "__data" / "libjupitermag" / "lib"
        for src_root in (data_root / "lib", data_root / "lib64", data_root / "bin"):
            if not src_root.is_dir():
                continue
            for src in src_root.rglob("*"):
                if not src.is_file() or src.suffix.lower() not in exts:
                    continue
                rel = src.relative_to(src_root)
                dst = dst_root / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


def getversion():
    """
    read the version string from __init__

    """
    # get the init file path
    thispath = os.path.abspath(os.path.dirname(__file__)) + "/"
    initfile = thispath + "JupiterMag/__init__.py"

    # read the file in
    f = open(initfile, "r", encoding="utf-8")
    lines = f.readlines()
    f.close()

    # search for the version
    version = "unknown"
    for line in lines:
        if "__version__" in line:
            s = line.split("=")
            version = s[-1].strip().strip('"').strip("'")
            break
    return version


version = getversion()

setup(
    name="JupiterMag",
    version=version,
    author="Matthew Knight James",
    author_email="mattkjames7@gmail.com",
    description="Some magnetic field models for Jupiter",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mattkjames7/JupiterMag",
    packages=find_packages(),
    package_data={
        "JupiterMag": [
            "__data/libjupitermag/lib/*.so",
            "__data/libjupitermag/lib64/*.so",
            "__data/libjupitermag/lib/*.dylib",
            "__data/libjupitermag/lib64/*.dylib",
            "__data/libjupitermag/lib/*.dll",
            "__data/libjupitermag/lib64/*.dll",
            "__data/libjupitermag/bin/*.dll",
        ]
    },
    cmdclass={"build_py": CustomBuild, "bdist_wheel": bdist_wheel},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: POSIX",
    ],
    install_requires=[
        "numpy",
        "matplotlib",
        "DateTimeTools",
        "RecarrayTools",
        "PyFileIO",
        "scipy",
    ],
    include_package_data=False,
    distclass=BinaryDistribution,
)

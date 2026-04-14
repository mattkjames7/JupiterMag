from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.dist import Distribution
import subprocess
import os
import platform
import shutil
from pathlib import Path

try:
    # Setuptools now provides bdist_wheel directly.
    from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
except ImportError:
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
        cmake_lists = lib_root / "CMakeLists.txt"

        if cmake_lists.is_file():
            self._build_with_cmake(lib_root)
        else:
            self._build_with_make(lib_root)

        # Fail fast if the native library was not produced.
        expected_lib = self._main_library_path(lib_root)
        if not expected_lib.is_file():
            raise RuntimeError(f"Native library was not built at {expected_lib}")

    def _build_with_cmake(self, lib_root: Path):
        build_cmd = self.get_finalized_command("build")
        build_dir = Path(build_cmd.build_temp) / "libjupitermag-cmake"
        build_dir.mkdir(parents=True, exist_ok=True)

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

    def _build_with_make(self, lib_root: Path):
        if platform.system() == "Windows":
            subprocess.check_call(["make", "-C", str(lib_root), "windows"], stderr=subprocess.STDOUT)
            return

        subprocess.check_call(["make", "-C", str(lib_root), "all"], stderr=subprocess.STDOUT)

    def _main_library_path(self, lib_root: Path) -> Path:
        ext = {
            "Windows": "dll",
            "Linux": "so",
            "Darwin": "dylib",
        }.get(platform.system())
        if ext is None:
            raise RuntimeError(f"Unsupported OS: {platform.system()}")
        return lib_root / "lib" / f"libjupitermag.{ext}"

    def copy_native_artifacts(self):
        """Copy generated native artifacts into build_lib for wheel packaging."""
        exts = {".dll", ".so", ".dylib"}
        root = Path(__file__).resolve().parent
        data_root = root / "JupiterMag" / "__data" / "libjupitermag"
        src_root = data_root / "lib"
        if not src_root.is_dir():
            return

        dst_root = Path(self.build_lib) / "JupiterMag" / "__data" / "libjupitermag" / "lib"
        for src in src_root.rglob("*"):
            if not src.is_file() or src.suffix.lower() not in exts:
                continue
            rel = src.relative_to(src_root)
            dst = dst_root / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)


setup(
    cmdclass={"build_py": CustomBuild, "bdist_wheel": bdist_wheel},
    distclass=BinaryDistribution,
)

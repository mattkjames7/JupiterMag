from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import subprocess
import os
import platform
import shutil
from pathlib import Path

class CustomBuild(build_py):
    def _log(self, message):
        self.announce(f"[JupiterMag build] {message}", level=2)

    def run(self):
        self._log("Starting CustomBuild.run")
        self.execute(self.target_build, ())
        build_py.run(self)
        self.copy_native_artifacts()
        self._log("Finished CustomBuild.run")

    def target_build(self):
        root = Path(__file__).resolve().parent
        lib_root = root / 'JupiterMag' / '__data' / 'libjupitermag'
        expected_lib = self._main_library_path(lib_root)

        self._log(f"Project root: {root}")
        self._log(f"Native lib root: {lib_root}")
        self._log(f"Expected main library: {expected_lib}")

        if platform.system() == 'Windows':
            cwd = os.getcwd()
            try:
                os.chdir(str(lib_root))
                self._log(f"Running Windows build script in: {os.getcwd()}")
                subprocess.check_call(['cmd', '/c', 'compile.bat'])
            finally:
                os.chdir(cwd)
        else:
            self._log("Running make for non-Windows native build")
            subprocess.check_call(['make', '-C', 'JupiterMag/__data/libjupitermag'])

        if (lib_root / 'lib').is_dir():
            self._log("Contents of native lib directory after build:")
            for p in sorted((lib_root / 'lib').rglob('*')):
                if p.is_file():
                    self._log(f"  {p}")

        # Fail fast if the native library was not produced.
        if not expected_lib.is_file():
            raise RuntimeError(
                f"Native library was not built at {expected_lib}"
            )

    def _main_library_path(self, lib_root: Path) -> Path:
        ext = {
            'Windows': 'dll',
            'Linux': 'so',
            'Darwin': 'dylib',
        }.get(platform.system())
        if ext is None:
            raise RuntimeError(f"Unsupported OS: {platform.system()}")
        return lib_root / 'lib' / f'libjupitermag.{ext}'

    def copy_native_artifacts(self):
        """Copy generated native artifacts into build_lib for wheel packaging."""
        root = Path(__file__).resolve().parent
        src_root = root / 'JupiterMag' / '__data' / 'libjupitermag' / 'lib'
        dst_root = Path(self.build_lib) / 'JupiterMag' / '__data' / 'libjupitermag' / 'lib'

        self._log(f"Copying native artifacts from {src_root} to {dst_root}")

        if not src_root.is_dir():
            self._log("Source native lib directory does not exist; skipping copy")
            return

        exts = {'.dll', '.so', '.dylib', '.a', '.lib'}
        for src in src_root.rglob('*'):
            if not src.is_file() or src.suffix.lower() not in exts:
                continue
            rel = src.relative_to(src_root)
            dst = dst_root / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            self._log(f"Copied artifact: {src} -> {dst}")


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

def getversion():
    '''
    read the version string from __init__
    
    '''
    #get the init file path
    thispath = os.path.abspath(os.path.dirname(__file__))+'/'
    initfile = thispath + 'JupiterMag/__init__.py'
    
    #read the file in
    f = open(initfile,'r',encoding='utf-8')
    lines = f.readlines()
    f.close()
    
    #search for the version
    version = 'unknown'
    for l in lines:
        if '__version__' in l:
            s = l.split('=')
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
    package_data={'JupiterMag': ['**/*']},
    cmdclass={'build_py': CustomBuild},  
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: POSIX",
    ],
    install_requires=[
        'numpy',
        'matplotlib',
        'DateTimeTools',
        'RecarrayTools',
        'PyFileIO',
        'scipy',
    ],
    include_package_data=True,
)




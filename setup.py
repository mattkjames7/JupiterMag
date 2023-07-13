from setuptools.command.install import install
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import subprocess
import os
import platform

class CustomBuild(build_py):
    def run(self):
        self.execute(self.target_build, ())
        build_py.run(self)

    def target_build(self):
        if platform.system() == 'Windows':
            cwd = os.getcwd()
            os.chdir('JupiterMag/__data/libjupitermag/')
            subprocess.check_call(['cmd','/c','compile.bat'])
            os.chdir(cwd)
        else:
            subprocess.check_call(['make', '-C', 'JupiterMag/__data/libjupitermag'])


with open("README.md", "r") as fh:
    long_description = fh.read()

def getversion():
	'''
	read the version string from __init__
	
	'''
	#get the init file path
	thispath = os.path.abspath(os.path.dirname(__file__))+'/'
	initfile = thispath + 'JupiterMag/__init__.py'
	
	#read the file in
	f = open(initfile,'r')
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




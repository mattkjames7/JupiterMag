# JupiterMag

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.6822191.svg)](https://doi.org/10.5281/zenodo.6822191)

Python wrapper for a collection of Jovian magnetic field models written in C++ (see [libjupitermag](https://github.com/mattkjames7/libjupitermag.git)).

This is part of a community code project : 

[Magnetospheres of the Outer Planets Group Community Code](https://lasp.colorado.edu/mop/missions/juno/community-code/)

**Journal Paper DOI**: [https://doi.org/10.1007/s11214-023-00961-3](https://doi.org/10.1007/s11214-023-00961-3)
(PDF via DOI, or [https://rdcu.be/c5I71](https://rdcu.be/c5I71), see [Journal Publication](README.md#journal-publication).)

**Authors**

- Matt James - University of Leicester

- Gabby Provan - University of Leicester

- Aneesah Kamran - University of Leicester

- Rob Wilson - LASP

- Marissa Vogt - Boston University

- Marty Brennan - NASA JPL

- Stan Cowley - University of Leicester

## 1 Requirements

For the Python code to run, the following Python packages would be required:

- NumPy

- Matplotlib

- DateTimeTools

- RecarrayTools

- PyFileIO

all of which would be installed automatically if using `pip`.

During installation, the C++ library which this module uses will be compiled.

### 1.1 Linux

JupiterMag was built and tested primarily using Linux Mint 20.3 (based on Ubuntu 20.04/Debian). To rebuild the code, ensure that `g++`, `make` and `ld` are installed.

### 1.2 Windows

This has been tested on Windows 10 (64-bit), other versions may also work. Requires `g++` and `ld` to work (these can be provided by TDM-GCC). This may or may not work with other compilers installed.

### 1.3 MacOS

This module has been tested on MacOS 11 Big Sur. It requires `g++`, `make` and `libtool` to recompile (provided by Xcode).

## 2 Installation

Install using `pip3`:

```bash
pip3 install JupiterMag --user
```

Download the latest release (on the right -> if you're viewing this on GitHub), then from within the directory where it was saved:

```bash
pip3 install JupiterMag-1.2.0.tar.gz --user
```

Or using this repo (replace "1.2.0" with the current version number):

```bash
#pull this repo
git clone https://github.com/mattkjames7/JupiterMag.git
cd JupiterMag

#update the submodule
git submodule update --init --recursive

#build the source distribution file
python3 setup.py sdist
#the output of the previous command should give some indication of 
#the current version number. If it's not obvious then do
# $ls dist/ to see what the latest version is
pip3 install dist/JupiterMag-1.2.0.tar.gz --user
```

I recommend installing `gcc` >= 9.3 (that's what this is tested with, earlier versions may not support the required features of C++).

This module should now work with both Windows and MacOS

### 2.1 Update an Existing Installation

To update an existing installation:

```bash
pip3 install JupiterMag --upgrade --user
```

Alternatively, uninstall then reinstall, e.g.:

```bash
pip3 uninstall JupiterMag
pip3 install JupiterMag --user
```

## 3 Usage

### 3.1 Internal Field

A number of internal field models are included (see [here](https://github.com/mattkjames7/libinternalfield/blob/main/README.md) for more information) and can be accessed via the ```JupiterMag.Internal``` submodule, e.g.:

```python
import JupiterMag as jm

#configure model to use VIP4 in polar coords (r,t,p)
jm.Internal.Config(Model="vip4",CartesianIn=False,CartesianOut=False)
Br,Bt,Bp = jm.Internal.Field(r,t,p)

#or use jrm33 in cartesian coordinates (x,y,z)
jm.Internal.Config(Model="jrm33",CartesianIn=True,CartesianOut=True)
Bx,By,Bz = jm.Internal.Field(x,y,z)
```

All coordinates are either in planetary radii (`x,y,z,r`) or radians (`t,p`). All Jovian models here use _R<sub>j</sub>_=71,492 km.

NOTE: figure 1 of the paper, which presents the radial components of the JRM33 model using degrees of 13 and 18 is corrupted, here is a clean version created using this code:

![](paper-fig1.png)

### 3.2 External Field

Currently the only external field source included is the Con2020 field (see [here](https://github.com/gabbyprovan/con2020.git) for the standalone Python code and [here](https://github.com/mattkjames7/libcon2020.git) for more information on the C++ code used here as part of libjupitermag), other models could be added in future.

This works in a similar way to the internal field, e.g.:

```python
#configure model
jm.Con2020.Config(equation_type='analytic')
Bx,By,Bz = jm.Con2020.Field(x,y,z)
```

### 3.3 Tracing

Field line tracing can be done using the `TraceField` object, e.g.

```python
import JupiterMag as jm

#configure external field model prior to tracing
#in this case using the analytic Con2020 model for speed
jm.Con2020.Config(equation_type='analytic')

#trace the field in both directions from a starting position
T = jm.TraceField(5.0,0.0,0.0,IntModel='jrm09',ExtModel='Con2020')
```

The above example will trace the field line from the Cartesian SIII position (5.0,0.0,0.0) (R<sub>j</sub>) in both directions until it reaches the planet using the JRM09 internal field model with the Con2020 external field model. The exact point at which the trace stops is either the surface of the planet or the radius L-MIC ionosphere, whichever has the shortest radial distance to the centre of Jupiter (see 3.3.2). The object  returned, `T`, is an instance of the `TraceField` class which contains the positions and magnetic field vectors at each step along the trace, along with some footprint coordinates and member functions which can be used for plotting.

#### 3.3.1 Plotting

The `TraceField` class has a few methods  which will create some simple plots of the traced field lines:
1. `PlotXZ()` - plot the field traces in the X-Z plane (SIII)
2. `PlotXY()` - plot the traces projected in the X-Y plane (SIII)
3. `PlotRhoZ()` - plot the field traces in ρ-Z plane, where ρ<sup>2</sup> = x<sup>2</sup> + y<sup>2</sup>

The example below can be used to compare field traces using just an internal field model (JRM33) with both internal and external field models (JRM33  + Con2020):

```python
import JupiterMag as jm
import numpy as np

#be sure to configure external field model prior to tracing
jm.Con2020.Config(equation_type='analytic')
#this may also become necessary with internal models in future, e.g.
#setting the model degree

#create some starting positions
n = 8
theta = (180.0 - np.linspace(22.5,35,n))*np.pi/180.0
r = np.ones(n)
x0 = r*np.sin(theta)
y0 = np.zeros(n)
z0 = r*np.cos(theta)

#create trace objects, pass starting position(s) x0,y0,z0
T0 = jm.TraceField(x0,y0,z0,Verbose=True,IntModel='jrm33',ExtModel='none')
T1 = jm.TraceField(x0,y0,z0,Verbose=True,IntModel='jrm33',ExtModel='Con2020')

#plot a trace
ax = T0.PlotRhoZ(label='JRM33',color='black')
ax = T1.PlotRhoZ(fig=ax,label='JRM33 + Con2020',color='red')

ax.set_xlim(-2.0,15.0)
ax.set_ylim(-6.0,6.0)
```

The resulting objects T0 and T1 store arrays of trace positions and magnetic field vectors along with a bunch of footprints.The above code produces a plot like this:

![](CompareTrace.png)


#### 3.3.2 Trace Footprints

The field traces each have an associated set of "footprints" which are stored
within the `TraceField` instance as a set of three `numpy.recarray`s called:

- `TraceField.surface` - locations where the field line intersects the
"surface" of Jupiter (oblate spheroid, equatorial radius = 71,492 km and polar 
radius = 66,854 km).
- `TraceField.ionosphere` - where the field lines map to the ionosphere as 
used by th L-MIC model (sphere of radius 67,354 km, 500 km above polar radius).
- `TraceField.equator` - location of the farthest point of the field trace from
the planet (unless the field line is open).

Every trace has an element of each of the above footprint arrays, but if the 
field is open, then some of the elements will be filled with `NAN`, indicating 
the absence of a footprint.

The `surface` and `ionosphere` arrays have the following fields:

- `xn3`, `yn3`, `zn3` - north footprint, SIII coords (R<sub>j</sub>)
- `xs3`, `ys3`, `zs3` - south footprint, SIII coords (R<sub>j</sub>)
- `xnm`, `ynm`, `znm` - north footprint, dipole coords (R<sub>j</sub>)
- `xsm`, `ysm`, `zsm` - south footprint, dipole coords (R<sub>j</sub>)
- `latn`, `lonn` - latitude and longitude, north footprint, SIII coords (°)
- `lats`, `lons` - latitude and longitude, south footprint, SIII coords (°)
- `mlatn`, `mlonn` - latitude and longitude, north footprint, dipole coords (°)
- `mlats`, `mlons` - latitude and longitude, south footprint, dipole coords (°)

The `equator` array has these fields:

- `x3`, `y3`, `z3` - SIII coords (R<sub>j</sub>)
- `xm`, `ym`, `zm` - dipole coords (R<sub>j</sub>)
- `mlone` - longitude in dipole coords (°)
- `fllen` - length of the field line (R<sub>j</sub>)
- `mshell` - distance from the centre of the planet to the farthest point along
the field trace (R<sub>j</sub>)



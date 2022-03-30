# JupiterMag

Jupiter magnetic field models with plenty of bugs.

## Installation

(not implemented yet) Install using `pip3`:

```bash
pip3 install JupiterMag --user
```

Or using this repo:

```bash
git clone https://github.com/mattkjames7/JupiterMag.git
cd JupiterMag
python3 setup.py bdist_wheel
pip3 install dist/JupiterMag-x.x.x-py3-none-any.whl --user
```

I recommend installing `gcc` >= 9.3 (that's what this is tested with, earlier versions may not support the required features of C++).

Windows and MacOS support coming soon...

## Usage

### Internal Field

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

All coordinates are either in planetary radii (`x,y,z,r`) or radians (`t,p`).

### External Field

Currently the only external field source included is the Con2020 field (see [here](https://github.com/mattkjames7/Con2020))

This works in a similar way to the internal field, e.g.:

```python
#configure model
jm.Con2020.Config(equation_type='analytic')
Bx,By,Bz = jm.Con2020.Field(x,y,z)
```

### Tracing

There is an object for field tracing `JupiterMag.TraceField`:

```python
import JupiterMag as jm
import numpy as np

#be sure to configure external field model prior to tracing
jm.Con2020.Config(equation_type='analytic')
#this may also become necessary with internal models in future, e.g.
#setting the model degree

#create some starting positions
n = 8
theta = (180.0 - np.linspace(21,35,n))*np.pi/180.0
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
	
ax.set_xlim(-2.0,25.0)
ax.set_ylim(-10.0,10.0)
```

The resulting objects `T0` and `T1` store arrays of trace positions and magnetic field vectors along with a bunch of footprints.

The above code produces a plot like this:

![](CompareTrace.png)

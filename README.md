# python-airsea

[![pip](http://img.shields.io/pypi/v/airsea.svg?style=flat)](https://pypi.python.org/pypi/airsea)
![pip-downloads](http://img.shields.io/pypi/dm/airsea.svg?style=flat)
[![build-status](https://github.com/pyoceans/python-airsea/actions/workflows/tests.yml/badge.svg?style=flat)](https://github.com/pyoceans/python-airsea/actions)
[![license](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/pyoceans/python-airsea/blob/master/LICENSE.txt)

A translation of the original AIRSEA-2.0 [MATLAB toolbox](https://github.com/sea-mat/air-sea) for calculating the properties of airsea fluxes.

## Installation

To get the latest version, install directly from this Github repository using pip
```
pip install git+https://github.com/pyoceans/python-airsea.git
```
or download the source and install the local version.
```
pip install .
```

Alternatively, install from pypi.
```
pip install airsea
```

## Examples

### Wind stress

Default Large and Pond (1981) formulation and 10m anemometer height
```python
wspd = [10., 0.2, 12., 20., 30., 50.] # wind speed (m/s)

from airsea import windstress as ws
tau = ws.stress(wspd)
```

Large and Pond (1981) formulation and 3m anemometer height
```python
tau = ws.stress(wspd, 3.)
```

Smith (1988) formulation with optional air temperature and density
```python
tau = ws.stress(wspd, drag='smith', Ta=9., rho_air=1.25)
```

### Atmosphere

Air density from temperature and relative humidity
```python
Ta = [5., 15., 23.1] # air temp (deg C)
rh = [95, 100, 50] # relative humidity (%)

from airsea import atmosphere as asea
rho_air = asea.air_dens(Ta, rh)
```

With optional pressure input
```python
Pa = 1000. # pressure (mb)

rho_air = asea.air_dens(Ta, rh, Pa)
```

## Other resources

The primary goal of this project is to provide functionality that is similar to commonly used Matlab functions. The science of air-sea bulk flux parameterization is continually evolving and updated algorithms are now available. Other related Python projects include
* [FluxEngine](https://github.com/oceanflux-ghg/FluxEngine)
* [pyCOARE](https://github.com/noaa-psd/pyCOARE)
* [xcoare](https://github.com/dcherian/xcoare)
* [AeroBulk](https://brodeau.github.io/aerobulk/)

## Contributing

Contributions are welcome. Feel free to raise an issue or submit a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests).
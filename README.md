# tri_dislocations_python

`tri_dislocations_python` (or `tde`) is a uncreatively-named port of Brendan
Meade's [tde][tde_bjm] MATLAB code for generating displacements, strains and
stresses from a triangular dislocation within an elastic halfspace. The
algorithms developed by Meade are written about in Meade (2007): "Algorithms
for the calculation of exact displacements, strains, and stresses for
triangular dislocation elements in a uniform elastic half space", *Computers &
Geosciences*.

The Python code here is only minimally modified from the MATLAB source. The
biggest changes are in making a central `tde` module with routines common
to both strain and displacement calculations, and moving the nitty gritty
strain and displacement functions to their own files, which are called by
`tde` routines.

Also, line breaks.

The code is not exhaustively tested but seems to match the results produced by
the MATLAB code (run in Octave) to the expected precision following ~50,000
sequential arithmetic operations on floats (i.e., still good enough for earth
science).


## Citations

First off, Brendan Meade did all of the work in developing the mathematical
model and writing the MATLAB code that this is based on.  Please cite Meade (2007)
as referenced in the first paragraph of this readme.

If you're really generous you can also cite this repository, using this DOI and
my name (Richard Styron).
[![DOI](https://zenodo.org/badge/41379836.svg)](https://zenodo.org/badge/latestdoi/41379836)

The citation would look something like:

Styron, R., 2019, tri_dislocations_python, *Zenodo*, https://github.com/cossatot/tri_dislocations_python, DOI: 10.5281/zenodo.3368513

## Usage

### Installation

`tde` only requires `numpy` and the standard Python library. It is tested on
Python 3.4 but should run on any system. It optionally uses `cython` which
greatly increases the speed of many of the strain calculations.

#### Building Cython files
run `python setup.py build_ext --inplace` in the outer `tri_dislocations_python`
directory.


### Installation

Run `python setup.py install` in the outer `tri_dislocations_python` directory
to install it system-wide. 

All the real code is in the `tde` folder. If system-wide installation is not
desired, copy that folder to the desired directory or temporarily add it to
the `$PYTHONPATH` or with `sys.path.append`.


### Usage

To calculate displacements, run
```python
import tde

U = tde.calc_tri_displacements(sx, sx, sz, x, y, z, pr, ss, ts, ds)
```
where `sx, sy, sz` are the 'station' coordinates or coordinates where the
calculations are to be made, `x, y, z` are arrays of the (x,y,z) coordinates
of the triangular dislocation's vertices, `pr` is the Poisson ratio, and
`ss, ts, ds` are the strike-slip displacement, tensile displacement, and
dip-slip displacement. `U` will be returned as a `dict` with `x`, `y`, and `z`
fields with those displacements for each station.

Similarly, strain is calculated as
```python
E = tde.calc_tri_strains(sx, sx, sz, x, y, z, pr, ss, ts, ds)
```
`E` here is a dict with each of the strain components.

Stress is calculated using `S` and material parameters:
```python
S = tde.strain_to_stress(E, lamda, mu)
```
`lamda` (NOT `lambda` which is a Python keyword) and `mu` are the Lame's
constants.


## In the Year 2000

The code is (likely) as functional as its source right now. Currently, a few
enhancements are planned:

- [ ] Optimization using [`numexpr`][numexpr] or [`numba`][numba], if those are
      installed on the user's system.

- [ ] Utilities for automatically running a collection of triangular
      displacements.


[tde_bjm]: https://github.com/brendanjmeade/tde
[numexpr]: https://github.com/pydata/numexpr
[numba]: http://numba.pydata.org/

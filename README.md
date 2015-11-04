turbinesFoam
============

[![Build Status](https://travis-ci.org/turbinesFoam/turbinesFoam.svg?branch=master)](https://travis-ci.org/turbinesFoam/turbinesFoam)
[![Stories in Ready](https://badge.waffle.io/turbinesfoam/turbinesfoam.png?label=ready&title=Ready)](https://waffle.io/turbinesfoam/turbinesfoam)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32471.svg)](http://dx.doi.org/10.5281/zenodo.32471)

`turbinesFoam` is a library for simulating wind and marine
hydrokinetic turbines in OpenFOAM (2.4.x) using the actuator line approach.

[![](https://cloud.githubusercontent.com/assets/4604869/10141523/f2e3ad9a-65da-11e5-971c-b736abd30c3b.png)](https://www.youtube.com/watch?v=THZvV4R1vow)


Status
------

This library is in development and is not yet fully functional.
See the [issue tracker](https://github.com/petebachant/turbinesFoam/issues)
for more details.
Pull requests are encouraged!

Also be sure to check out the
[development snapshot videos on YouTube](https://www.youtube.com/playlist?list=PLOlLyh5gytG8n8D3V1lDeZ3e9fJf9ux-e).


Features
--------

`fvOptions` classes for adding actuator lines and turbines constructed from
actuator lines to any compatible solver or turbulence model, e.g.,
`simpleFoam`, `pimpleFoam`, `interFoam`).


Installation
------------

```bash
cd $WM_PROJECT_USER_DIR
git clone https://github.com/turbinesFoam/turbinesFoam.git
cd turbinesFoam
./Allwmake
```


Usage
-----
There are tutorials located in `turbinesFoam/tutorials`


How to cite
-----------
```bibtex
@misc{Bachant2015-turbinesFoam,
  author       = {Peter Bachant and Martin Wosnik},
  title        = {turbinesFoam: v0.0.4},
  month        = October,
  year         = 2015,
  doi          = {10.5281/zenodo.32471},
  url          = {http://dx.doi.org/10.5281/zenodo.32471},
  howpublished = {Zenodo. \url{http://dx.doi.org/10.5281/zenodo.32471}}
}
```


License
-------

See LICENSE.

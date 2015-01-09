turbinesFoam
============

`turbinesFoam` is a library for parameterizing wind and marine 
hydrokinetic turbines in OpenFOAM. Some components for solvers 
and turbine models have been taken from NREL's SOWFA and Offwind.

Status
------

This library is in heavy development and is not yet fully functional.
See the [issue tracker](https://github.com/petebachant/turbinesFoam/issues)
for more details. 
Pull requests are encouraged!

Features
--------
`fvOptions` classes for adding turbine effects to any solver or turbulence model
that accepts these (e.g. `simpleFoam`, `pimpleFoam`, `interFoam`). 


Compiling
---------

```
cd $WM_PROJECT_USER_DIR
git clone https://github.com/petebachant/turbinesFoam.git
cd turbinesFoam
./Allwmake
```

Usage
-----
There are tutorials located in `turbinesFoam/tutorials`

License
-------

See LICENSE.

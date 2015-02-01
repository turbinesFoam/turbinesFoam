turbinesFoam
============

`turbinesFoam` is a library for simulating wind and marine 
hydrokinetic turbines in OpenFOAM (2.3.x). Some components for solvers 
and turbine models have been taken from NREL's SOWFA and Offwind.

Status
------

This library is in heavy development and is not yet fully functional.
See the [issue tracker](https://github.com/petebachant/turbinesFoam/issues)
for more details. 
Pull requests are encouraged!

Also be sure to check out the 
[development snapshot videos on YouTube](https://www.youtube.com/playlist?list=PLOlLyh5gytG8n8D3V1lDeZ3e9fJf9ux-e).

Features
--------
`fvOptions` classes for adding turbine effects (actuator line model)
to any solver or turbulence model
that accepts these (e.g. `simpleFoam`, `pimpleFoam`, `interFoam`). 


Compiling
---------

```
cd $WM_PROJECT_USER_DIR
git clone https://github.com/turbinesFoam/turbinesFoam.git
cd turbinesFoam
./Allwmake
```

Usage
-----
There are tutorials located in `turbinesFoam/tutorials`

License
-------

See LICENSE.

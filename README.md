turbinesFoam
============

Right now, this library is just a dream. Don't expect anything to work and please
consider helping with the project as I am a total beginner with C++! 

Someday, this will hopefully be a library for parameterizing wind and marine 
hydrokinetic turbines in OpenFOAM. So far, some components for solvers 
and turbine models have been taken from NREL's SOWFA and Offwind.

To-do
-----
  - [ ] Write axial-flow turbine actuator line model as `fvOption`, so custom
        solvers are not necessary.
  - [ ] Write cross-flow turbine actuator line model as `fvOption`.
  - [ ] Write cross-flow turbine actuator surface model as `fvOption`.


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

turbinesFoam
============

A library for actuator line turbine simulations in OpenFOAM.
Some components were taken from SOWFA and Offwind.


To-do
-----
  - [ ] Use fvOptions for turbine models, so solvers do not need to be modified?
  - [ ] Write cross-flow turbine actuator line model
  - [ ] Use turbines as sources in turbulence model equations?


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

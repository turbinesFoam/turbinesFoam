turbinesFoam
============

A library for wind and marine hydrokinetic turbine parameterized simulations 
in OpenFOAM. Some components for solvers and turbine models were taken from 
SOWFA and Offwind.

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

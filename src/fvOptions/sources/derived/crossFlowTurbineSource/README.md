crossFlowTurbineSource
======================

Modified `actuationDiskSource` to better represent cross-flow turbines. 

Compilation/installation
------------------------
Execute `wmake libso`.

Usage
-----

In `system/controlDict`:

```
libs
(
    "libOpenFOAM.so"
    "libincompressibleTurbulenceModel.so"
    "libincompressibleRASModels.so"
    "libmyFvOptions.so"
);
```



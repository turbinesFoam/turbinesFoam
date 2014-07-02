crossFlowTurbineActuatorLineSource
==================================

Modified `rotorDiskSource` to represent cross-flow turbines. 

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

Design
------
* Create a cellZone in the mesh that will locate the turbine.
* Blades will sweep along the outer surface where the surface is not parallel
  to the axis of rotation. This means any frontal area can be used. 
* What if surface is not symmetrical about axis, e.g., a cube? Blade elements
  could always move tangential to surface.
* Blade element input data:

|   `profile`  | `axialDistance` | `radius` | `chord` | `chordMount` | `pitch` |
|--------------|-----------------|----------|---------|--------------|---------|
| `"naca0020"` |       0.0       |    0.5   |   0.14  |     0.5      |   0.0   |

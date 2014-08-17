crossFlowTurbineALSource
========================

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
* Sample input dict:

```
crossFlowTurbineALCoeffs
{
    fieldNames	    (U);
    nBlades	        3;
    tipEffect	    1.0; 

    inletFlowType   local;  // inlet flow type specification
    inletVelocity   (0 1 0);
    upstreamPoint   (-1.0, 0, 0); // Absolute or relative?
    tipSpeedRatio   1.9;

    geometryMode    specified;   // geometry specification (???)
    origin          (0 0 0);
    axis            (0 0 1);

    blades
    {
        blade0
        {
            elementData
            ( // axialDistance, radius, azimuth, profile, chord, chordMount, pitch
                (0.0 0.5 0.0 NACA0020 0.14 0.5 0.0)
                (1.0 0.5 0.0 NACA0020 0.14 0.5 0.0)
            );
        }
        blade1
        {
            elementData
            ( // axialDistance, radius, azimuth, profile, chord, chordMount, pitch
                (0.0 0.5 120.0 NACA0020 0.14 0.5 0.0)
                (1.0 0.5 120.0 NACA0020 0.14 0.5 0.0)
            );
        }
        blade1
        {
            elementData
            ( // axialDistance, radius, azimuth, profile, chord, chordMount, pitch
                (0.0 0.5 240.0 NACA0020 0.14 0.5 0.0)
                (1.0 0.5 240.0 NACA0020 0.14 0.5 0.0)
            );
        };
    };
    
    profiles
    {
        NACA0020
        {
	    type 	lookup;
	    data 	(
			(-90 0.21 1.45)
			(-18 0.21 1.45)
			(-16 0.165 1.3)
			(-14 0.125 1.1)
			(-12 0.092 0.95)
			(-10 0.07 0.8)
			(-8 0.05 0.64)
			(-6 0.04 0.5)
			(-4 0.028 0.32)
			(-2 0.022 0.18)
			(0 0.02 0)
			(2 0.022 0.18)
			(4 0.028 0.32)
			(6 0.04 0.5)
			(8 0.05 0.64)
			(10 0.07 0.8)
			(12 0.092 0.95)
			(14 0.125 1.1)
			(16 0.165 1.3)
			(18 0.21 1.45)
			(90 0.21 1.45)
			);
        };
    };
}
```

  * `axialDistance` is defined in meters as the linear distance from the
    most negative point of the axis.
  * 

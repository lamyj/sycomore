Common features
===============

.. _units:

Units
-----

MRI simulations deal with various quantities: times, frequencies and angular frequencies, magnetic field strength, gradient moments, etc. All those quantities use different usual prefixes, but not always the same: relaxation time and sequence timing are usually expressed in milliseconds, while the gyromagnetic ratio is expressed in *MHz/T* (i.e. inverse of microseconds) gradients moments are expressed in *mT/m⋅ms* and slice thickness in millimeters (on clinical scanners) or micrometers (on pre-clinical scanners). This wealth of units makes it very easy to get quantities wrong by a factor of 1000 or more.

Sycomore provides a unit system so that users do not have to convert their quantities to a specific unit. Units may be declared by multiplying or dividing by the unit name (e.g. ``500*ms``). Those two syntaxes can be mixed in order to use more complex units (e.g. ``267.522*MHz/T``). Unit objects follow the usual arithmetic rules, and all SI `base units`_, `derived units`_ and `prefixes`_ are defined. 

Quantity (:cpp:class:`C++ <sycomore::Quantity>`, :attr:`Python <sycomore.Quantity>`) objects contain their value in the base SI unit and may be converted to a compatible unit. Common arithmetic operations (addition, subtraction, multiplication, division, power) are implemented, as well as most of `numpy numeric functions`_. The following code sample summarizes these features.

.. note:: The *micro* prefix is *u*, not *μ*, in order to keep ASCII names


.. code:: python

    import sycomore
    from sycomore.units import *
    
    # Combination of base units
    diffusion_coefficient = 0.89*um**2/ms
    
    # Magnitude of the quantity, in SI unit
    length = 180*cm
    length_in_meters = length.magnitude # Equals to 1.8
    
    # Conversion: magnitude of the quantity, in specified unit
    duration = 1*h
    duration_in_seconds = duration.convert_to(s) # Equals to 3600

.. note:: Take caution when importing all units name (``from sycomore.units import *``), as some names may clash with your code. In a long module with many variables, it is better to import only the required units (e.g. ``from sycomore.units import mT, m, ms``)

Species
-------

A species is characterized by its relaxation rates (|R1|, |R2| and |R2'|), its diffusivity *D* and its relative resonance frequency Δω. It is described by the Species (:cpp:class:`C++ <sycomore::Species>`, :attr:`Python <sycomore.Species>`) class. |R1| and |R2| are mandatory parameters, all others are optional and are equal to 0 if unspecified. Using the :ref:`unit system<units>`, relaxation rates and relaxation times may be used interchangeably.

.. code:: python

    import sycomore
    from sycomore.units import *
    
    # Create a Species from either relaxation times, relaxation rates or both
    species = sycomore.Species(1000*ms, 100*ms)
    species = sycomore.Species(1*Hz, 10*Hz)
    species = sycomore.Species(1000*ms, 10*Hz)

The diffusivity can be assigned either as a scalar (for isotropic diffusion) or as a tensor (for anisotropic diffusion), but will always be returned as a tensor:

.. code:: python

    import sycomore
    from sycomore.units import *
    
    species = sycomore.Species(1000*ms, 100*ms)
    # Assign the diffusion coefficient as a scalar
    species.D = 3*um**2/s
    # The diffusion coefficient is stored on the diagonal of the tensor
    print(species.D[0])
    
    # Assign the diffusion coefficient as a tensor
    species.D = [
        3*um**2/s, 0*um**2/s, 0*um**2/s,
        0*um**2/s, 2*um**2/s, 0*um**2/s,
        0*um**2/s, 0*um**2/s, 1*um**2/s]
    print(species.D)

.. code::

    3e-12 [ L^2 T^-1 ]
    (3e-12 [ L^2 T^-1 ] 0 [ L^2 T^-1 ] 0 [ L^2 T^-1 ] 0 [ L^2 T^-1 ] 2e-12
    [ L^2 T^-1 ] 0 [ L^2 T^-1 ] 0 [ L^2 T^-1 ] 0 [ L^2 T^-1 ] 1e-12 [ L^2
    T^-1 ])

Time intervals
--------------

A time interval (:cpp:class:`C++ <sycomore::TimeInterval>`, :attr:`Python <sycomore.TimeInterval>`) is specified by its duration and an optional magnetic field gradient. The gradient can be either as a scalar or as a 3D array, and can describe the amplitude (in *T/m*), the area (in *T/m\*s*) or the dephasing (in *rad/m*). A :func:`sycomore.TimeInterval.set_gradient` function is available for generic modification of the gradient.

.. code:: python

    import sycomore
    from sycomore.units import *
    
    # Scalar gradient, defined by its amplitude
    interval = sycomore.TimeInterval(1*ms, 20*mT/m)
    print(
        interval.duration, "\n",
        interval.gradient_amplitude, "\n",
        interval.gradient_area/interval.duration, "\n",
        interval.gradient_dephasing/(sycomore.gamma*interval.duration))

.. code::

    0.001 [ T ]
     (0.02 [ L^-1 M T^-2 I^-1 ] 0.02 [ L^-1 M T^-2 I^-1 ] 0.02 [ L^-1 M
    T^-2 I^-1 ])
     (0.02 [ L^-1 M T^-2 I^-1 ] 0.02 [ L^-1 M T^-2 I^-1 ] 0.02 [ L^-1 M
    T^-2 I^-1 ])
     (0.02 [ L^-1 M T^-2 I^-1 ] 0.02 [ L^-1 M T^-2 I^-1 ] 0.02 [ L^-1 M
    T^-2 I^-1 ])

.. |R1| replace:: R\ :sub:`1`
.. |R2| replace:: R\ :sub:`2`
.. |R2'| replace:: R\ :sub:`2`:sup:`'`

.. _base units: https://en.wikipedia.org/wiki/SI_base_unit
.. _derived units: https://en.wikipedia.org/wiki/SI_derived_unit
.. _numpy numeric functions: https://docs.scipy.org/doc/numpy/reference/ufuncs.html#available-ufuncs
.. _prefixes: https://en.wikipedia.org/wiki/Metric_prefix

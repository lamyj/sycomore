Extended Phase Graphs (EPG)
===========================

Extended Phase Graphs (EPG) are first described in Hennig's 1988 paper `Multiecho imaging sequences with low refocusing flip angles`_ -- although the term "Extended Phase Graph" appears in a later paper (`Echoes—how to generate, recognize, use or avoid them`_) -- is an invaluable tool for fast simulation of MRI sequences. It allows simulation of all main phenomena, including `diffusion`_, `motion`_, and `magnetization transfer`_.

The EPG model is based on the spatial Fourier transform of the isochromats present in a voxel. Originally developed for sequences with a very regular time course (e.g. RARE or bSSFP), it has been later extended for `non-regular sequences`_

Sycomore provides three EPG models: 

- a one-dimensional regular model where all gradient moments throughout the simulation are assumed equal, 
- a one-dimensional discrete model where gradients moments may differ, their values being discretized (i.e. binned) to avoid numerical instabilities along the simulation,
- a three-dimensional discrete model, where gradients can be specified in three dimensions.

All three models share a common API, and can handle single-pool or two-pools (exchange and MT) variants.

To create a single-pool model, a single species is sufficient. Additionally, the initial magnetization :math:`M_0` can be specified. If missing, it defaults to :math:`[0,0,1]`.

.. code:: python
   
   species_a = sycomore.Species(779*ms, 45*ms)
   single_pool_model = sycomore.epg.Discrete(species_a)

For an exchange model two species, their initial magnetizations, and the exchange rate from pool *a* to pool *b* must be provided. An optional frequency offset of pool *b* relative to pool *a* may be provided, default to 0 Hz.

.. code:: python
   
   species_b = sycomore.Species(100*ms, 20*ms)
   M0_a = sycomore.Array[float](0, 0, 0.8)
   M0_b = sycomore.Array[float](0, 0, 0.2)
   k_a = 2*Hz
   exchange_model = sycomore.epg.Discrete(species_a, species_b, M0_a, M0_b, k_a)

In an MT model, pool *b* is assumed to have a very short :math:`T_2` and thus no transversal magnetization: one species, the :math:`R_1` or :math:`T_1` of pool *b*, both initial magnetizations, and the exchange rate from pool *a* to pool *b* must be provided.

.. code:: python
   
   R1_b = 100*ms
   mt_model = sycomore.epg.Discrete(species_a, R1_b, M0_a, M0_b, k_a)

Once the model is created, the main features are:

- Simulation of an RF pulse using the function :code:`apply_pulse`. For a single pool model, the flip angle must be provided, and the phase is optional, defaulting to 0°. For an exchange model, a second pair of flip angle and phase may be provided -- if missing, the same pulse will be applied to the two pools. For an MT model, the phase and direct saturation are mandatory.
- Simulation of a time interval (relaxation, diffusion, exchange, etc.) using the function :code:`apply_time_interval`. Refer to the model-specific pages for details.
- Accessing the non-dephased magnetization using the :code:`echo` property
- Accessing the orders and the states using the :code:`orders` and :code:`states` properties.

.. toctree::
   :maxdepth: 1
   :caption: Contents:
   
   regular.rst
   discrete.rst
   discrete_3d.rst

.. _diffusion: 
.. _non-regular sequences: https://doi.org/10.1016/j.jmr.2010.05.011
.. _Echoes—how to generate, recognize, use or avoid them: https://doi.org/10.1002/cmr.1820030302
.. _magnetization transfer: https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27040
.. _motion: https://doi.org/10.1002/jmri.24619
.. _Multiecho imaging sequences with low refocusing flip angles: https://doi.org/10.1016/0022-2364(88)90128-X

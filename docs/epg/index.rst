Extended Phase Graphs (EPG)
===========================

Extended Phase Graphs (EPG) are first described in Hennig's 1988 paper `Multiecho imaging sequences with low refocusing flip angles`_ -- although the term "Extended Phase Graph" appears in a later paper (`Echoes—how to generate, recognize, use or avoid them`_) -- is an invaluable tool for fast simulation of MRI sequences. It allows simulation of all main phenomena, including `diffusion`_, `motion`_, and `magnetization transfer`_.

The EPG model is based on the spatial Fourier transform of the isochromats present in a voxel. Originally developped for sequences with a very regular time course (e.g. RARE or bSSFP), it has been later extended for `non-regular sequences`_

Sycomore provides three EPG models: 

- a one-dimensional regular model where all gradient moments throughout the simulation are assumed equal, 
- a one-dimensional discrete model where gradients moments may differ, their values being discretized (i.e. binned) to avoid numerical instabilities along the simulation,
- a three-dimensional discrete model, where gradients can be specified in three dimensions.

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

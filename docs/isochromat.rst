**todo**
use `..tab` as in https://raw.githubusercontent.com/pypa/setuptools/main/docs/userguide/quickstart.rst to display both Python and C++ code

separate userguide from api reference

do two references: python and c++ (https://breathe.readthedocs.io/en/latest/ ?)

Isochromat simulation
=====================

The isochromat simulation in Sycomore is based on Hargreaves's 2001 paper `Characterization and Reduction of the Transient Response in Steady-State MR Imaging`_ where limiting cases of the full spin behavior (instantaneous RF pulse, "pure" relaxation, "pure" precession, etc.) are expressed as matrices. The matrices representing building blocks of a sequence are multiplied amongst them, forming a single matrix operator for a repetition. Iterating this process yields a fast simulation of the evolution of a single isochromat and the eigenanalysis of the resulting matrix yields important insights on the steady-state of the sequence.

Homogeneous coordinates
-----------------------

An `homogeneous form`_ of all matrices and vectors is used to obtain a purely multiplicative form of the matrix operators. Using homogenous coordinates (also called projective coordinates), a vector in :math:`\mathbb{R}^3` with Cartesian coordinates :math:`(x_c, y_c, z_c)` is represented as a *set* of 4D vectors in :math:`\mathbb{PR}^3`, the three-dimensional projective space, :math:`(x_p, y_p, z_p, w)` so that :math:`(x_c, y_c, z_c) = (x_p/w, y_p/w, z_p/w)`. In this representation, any geometric transformation which can be expressed as a 3×3 matrix :math:`M` can also be represented as a 4×4 homogeneous matrix:

.. math::
  
  \begin{bmatrix}
    M & \begin{matrix} 0 \\ 0 \\ 0 \end{matrix} \\
    \begin{matrix} 0 & 0 & 0 \end{matrix} & 1
  \end{bmatrix}

Moreover, a translation :math:`\mathbf{T} = (T_x, T_y, T_z)` can also be represented as a 4×4 homogeneous matrix:

.. math::
  
  \begin{bmatrix}
    1 & 0 & 0 & T_x \\
    0 & 1 & 0 & T_y \\
    0 & 0 & 1 & T_z \\
    0 & 0 & 0 & 1
  \end{bmatrix}

The composition of the geometric transforms described by those homogeneous matrices is given by multiplying the matrices, as is the case in Euclidean coordinates.

Usage
-----

The basic operators of Bloch simulation are :func:`sycomore.bloch.pulse` and :func:`sycomore.bloch.time_interval`. They both use the :class:`sycomore.Species` class and the :ref:`units<units>` system. The following code sample shows the simulation of a saturation-recuperation experiment: the system is idle for the first 100 ms, a 60° pulse is applied, and the system then relaxes for 1 s. The time step is 10 ms.


.. code:: python

    import matplotlib.pyplot
    import numpy
    import sycomore
    from sycomore.units import *
    
    species = sycomore.Species(1000*ms, 100*ms)
    flip_angle = 60*deg
    
    idle = sycomore.bloch.time_interval(species, 10*ms)
    pulse = sycomore.bloch.pulse(flip_angle)
    
    t = 0*s
    M = numpy.array([0,0,1,1])
    
    record = [[t, M[:3]/M[3]]]
    for _ in range(10):
        t = t+10*ms
        M = idle @ M
        record.append([t, M[:3]/M[3]])
    
    M = pulse @ M
    record.append([t, M[:3]/M[3]])
    
    for _ in range(100):
        t = t+10*ms
        M = idle @ M
        record.append([t, M[:3]/M[3]])
    
    time, magnetization = list(zip(*record))
    magnetization = numpy.array(magnetization)
    
    x_axis = [x.convert_to(ms) for x in time]
    matplotlib.pyplot.plot(
        x_axis, numpy.linalg.norm(magnetization[:, :2], axis=-1), label="$M_\perp$")
    matplotlib.pyplot.plot(x_axis, magnetization[:, 2], label="$M_z$")
    matplotlib.pyplot.xlim(0)
    matplotlib.pyplot.ylim(-0.02)
    matplotlib.pyplot.xlabel("Time (ms)")
    matplotlib.pyplot.ylabel("$M/M_0$")
    matplotlib.pyplot.legend()
    matplotlib.pyplot.tight_layout()


.. figure:: ./bloch_saturation_recuperation_1.png
   :width: 15 cm

   Saturation-recuperation using Bloch simulation



Reference
---------


.. _Characterization and Reduction of the Transient Response in Steady-State MR Imaging: https://doi.org/10.1002/mrm.1170
.. _homogeneous form: https://en.wikipedia.org/wiki/Homogeneous_coordinates

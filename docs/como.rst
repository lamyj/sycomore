Configuration Model
===================

The Configuration Model, described in `Ganter's 2018 communication`_ is a microscopic alternative to the macroscopic EPG model in which each different time interval is represented as a dimension in the model. Compared to EPG, the Configuration Model provides easier quantification of susceptibility effects.

The following code sample simulates a GRE sequence in which the readout is located at the middle of the TR.

.. code-block:: python

  import math
  import sycomore
  from sycomore.units import *

  species = sycomore.Species(1000*ms, 100*ms, 0.89*um*um/ms)
  m0 = sycomore.Magnetization(0, 0, 1)

  flip_angle = 40*deg
  TR = 500*ms

  TR_count = 10

  model = sycomore.como.Model(
      species, m0,
      [["half_echo", sycomore.TimeInterval(TR/2.)]])

  magnetization = []
  for i in range(TR_count):
      model.apply_pulse(
          sycomore.Pulse(flip_angle, (math.pi/3+(i%2)*math.pi)*rad))
      model.apply_time_interval("half_echo")
      magnetization.append(model.isochromat())
      model.apply_time_interval("half_echo")

Reference
---------

.. autoclass:: sycomore.como.Model
  :members:

.. _Ganter's 2018 communication: http://archive.ismrm.org/2018/5663.html

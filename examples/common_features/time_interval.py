import sycomore
from sycomore.units import *

# Scalar gradient, defined by its amplitude
interval = sycomore.TimeInterval(1*ms, 20*mT/m)
print(
    interval.duration,
    interval.gradient_amplitude,
    interval.gradient_area/interval.duration,
    interval.gradient_dephasing/(sycomore.gamma*interval.duration),
    sep="\n")

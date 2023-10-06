import numpy
import sycomore
from sycomore.units import *

species = sycomore.Species(1000*ms, 100*ms)
TE = 4*ms
train_length = 40

model = sycomore.epg.Regular(species)
signal = numpy.zeros(train_length, dtype=complex)

model.apply_pulse(90*deg)
for echo in range(train_length):
    model.apply_time_interval(TE/2)
    model.apply_pulse(180*deg)
    model.apply_time_interval(TE/2)
    
    signal[echo] = model.echo

import sys
import matplotlib.pyplot

destination, = sys.argv[1:]

times = [(1+x)*TE for x in range(train_length)]
x_axis = [x.convert_to(ms) for x in times]

figure, plot = matplotlib.pyplot.subplots(tight_layout=True, figsize=(8, 6))
magnitude = numpy.abs(signal)
plot.plot(x_axis, magnitude, ".", label="Simulated")
plot.plot(
    x_axis, [numpy.exp(-(x*species.R2).magnitude) for x in times],
    label="$T_2$ decay")
plot.set(xlabel = "Time (ms)", ylabel = r"$M_\perp/M_0$ (unitless)")
plot.legend()
figure.savefig(destination)

import sys

import sycomore
from sycomore.units import *

destination, = sys.argv[1:]

M0 = [0,0,1]
positions = [[0*mm, 0*mm, 0*mm]]
T1, T2 = 1000*ms, 100*ms
flip_angle = 60*deg
step = 10*ms

model = sycomore.isochromat.Model(T1, T2, M0, positions)

pulse = model.build_pulse(flip_angle)
idle = model.build_time_interval(10*ms)

record = [[0*s, model.magnetization[0]]]
for _ in range(10):
    model.apply(idle)
    record.append([record[-1][0]+step, model.magnetization[0]])

model.apply(pulse)
record.append([record[-1][0], model.magnetization[0]])

for _ in range(100):
    model.apply(idle)
    record.append([record[-1][0]+step, model.magnetization[0]])

import matplotlib.pyplot
import numpy

time, magnetization = list(zip(*record))
magnetization = numpy.array(magnetization)

x_axis = [x.convert_to(ms) for x in time]
figure, plot = matplotlib.pyplot.subplots(tight_layout=True, figsize=(8, 6))
plot.plot(
    x_axis, numpy.linalg.norm(magnetization[:, :2], axis=-1), label="$M_\perp$")
plot.plot(x_axis, magnetization[:, 2], label="$M_z$")
plot.set(xlim=0, ylim=-0.02, xlabel="Time (ms)", ylabel="$M/M_0$")
matplotlib.pyplot.legend()
matplotlib.pyplot.tight_layout()

figure.savefig(destination)

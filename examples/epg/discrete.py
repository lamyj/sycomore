import sys
import time

import numpy
import sycomore
from sycomore.units import *

destination, = sys.argv[1:]

species = sycomore.Species(1000*ms, 100*ms, 1700*um**2/s)

# Sequence parameters
flip_angle = 20*deg
TE = 2*ms, 8*ms
TR = 10*ms
slice_thickness = 1*mm
tau_readout = 1*ms
q = 70/cm

# Motion to k-space extremity and its associated gradient amplitude
k_max = 0.5 * 2*numpy.pi/slice_thickness
G_readout = k_max / sycomore.gamma / (tau_readout/2)

# Diffusion gradient
tau_diffusion = (TE[1]-tau_readout/2) - (TE[0]+tau_readout/2)
G_diffusion = q/(sycomore.gamma_bar*tau_diffusion)

models = [sycomore.epg.Discrete(species), sycomore.epg.Discrete(species)]
models[1].threshold = 1e-6

repetitions = int((4*species.T1/TR))

S_plus = numpy.zeros((len(models), repetitions), dtype=complex)
S_minus = numpy.zeros((len(models), repetitions), dtype=complex)
for index, model in enumerate(models):
    begin = time.time()
    
    for repetition in range(repetitions):
        # RF-pulse and idle until the first read-out 
        model.apply_pulse(flip_angle)
        model.apply_time_interval(TE[0]-tau_readout)
        
        # First echo
        model.apply_time_interval(tau_readout/2, -G_readout)
        model.apply_time_interval(tau_readout/2, G_readout)
        S_plus[index, repetition] = model.echo
        model.apply_time_interval(tau_readout/2, G_readout)
        
        # Diffusion gradient between the two echoes
        model.apply_time_interval(tau_diffusion, G_diffusion)
        
        # Second echo
        model.apply_time_interval(tau_readout/2, G_readout)
        S_minus[index, repetition] = model.echo
        model.apply_time_interval(tau_readout/2, G_readout)
        model.apply_time_interval(tau_readout/2, -G_readout)
        
        # Idle until the end of the TR
        model.apply_time_interval(TR-TE[1]-tau_readout)
        
        # Make sure the sequence timing is correct
        if repetition == 0:
            assert((model.elapsed-TR)/TR < 1e-6)
    
    end = time.time()
    print(
        "Threshold:", model.threshold, len(model), "orders",
        1e3*(end-begin), "ms")

import matplotlib.pyplot

figure, plot = matplotlib.pyplot.subplots(
    1, 1, sharex=True, tight_layout=True, figsize=(6, 6))

plot.plot(range(repetitions), numpy.abs(S_plus[0]), "-C0", label="$S^+$")
plot.plot(range(repetitions), numpy.abs(S_minus[0]), "-C1", label="$S^-$")
plot.plot(range(repetitions)[::15], numpy.abs(S_plus[1, ::15]), "xC0")
plot.plot(range(repetitions)[::15], numpy.abs(S_minus[1, ::15]), "xC1")

plot.plot([], [], "-k", label="No threshold")
plot.plot([], [], "xk", label=f"Threshold={models[1].threshold}")

plot.set(
    xlim=(0, repetitions-1), ylim=0,
    xlabel="Repetition", ylabel="Magnitude (a.u.)")
# for index, (s, phi) in enumerate(zip(signals, phase_steps)):
#     plots[0].plot(
#         range(repetitions), numpy.abs(s), f"C{index}", lw=1,
#         label=rf"$\phi_{{inc}} = {phi.convert_to(deg):.0f}°$")
#     plots[1].plot(range(repetitions), numpy.angle(s), f"C{index}", lw=1)
# 
# plots[0].set(ylim=0, ylabel="Magnitude (a.u.)")
# plots[1].set(
#     xlim=(0, repetitions-1), ylim=(-numpy.pi, +numpy.pi),
#     xlabel="Repetition", ylabel="Phase (rad)")
plot.legend()

figure.savefig(destination)

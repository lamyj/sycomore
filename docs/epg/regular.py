import sys

import numpy
import sycomore
from sycomore.units import *

destination, = sys.argv[1:]

species = sycomore.Species(1000*ms, 1000*ms)

# Sequence parameters
flip_angle=30*deg
TE = 5*ms
TR = 25*ms
phase_steps = [0*deg, 90*deg, 117*deg, 180*deg]
slice_thickness = 1*mm
tau_readout = 1*ms
repetitions = int((5*species.T1/TR))

# Motion to k-space extremity and its associated gradient amplitude
k_max = 0.5 * 2*numpy.pi/slice_thickness
G = k_max / sycomore.gamma / (tau_readout/2)

models = [
    sycomore.epg.Regular(species, unit_dephasing=k_max) for _ in phase_steps]

signals = numpy.zeros((len(models), repetitions), dtype=complex)
for r in range(0, repetitions):
    for index, (phase_step, model) in enumerate(zip(phase_steps, models)):
        phase = (phase_step * 1/2*(r+1)*r)
        
        # RF-pulse and idle until the readout
        model.apply_pulse(flip_angle, phase)
        model.apply_time_interval(TE-tau_readout)
        
        # Readout prephasing and first half of the readout
        model.apply_time_interval(-G, tau_readout/2)
        model.apply_time_interval(+G, tau_readout/2)
        
        # Echo at the center of the readout, cancel the phase imparted by the
        # RF-spoiling
        signals[index, r] = model.echo * numpy.exp(-1j*phase.convert_to(rad))
        
        # Second half of the readout, idle until the end of the TR
        model.apply_time_interval(+G, tau_readout/2)
        model.apply_time_interval(TR-TE-tau_readout/2)

import matplotlib.pyplot

figure, plots = matplotlib.pyplot.subplots(
    2, 1, sharex=True, tight_layout=True, figsize=(8, 8))
for index, (s, phi) in enumerate(zip(signals, phase_steps)):
    plots[0].plot(
        range(repetitions), numpy.abs(s), f"C{index}",
        label=rf"$\phi_{{inc}} = {phi.convert_to(deg):.0f}°$")
    plots[1].plot(range(repetitions), numpy.angle(s), f"C{index}")

plots[0].set(ylim=0, ylabel="Magnitude (a.u.)")
plots[1].set(
    xlim=(0, repetitions-1), ylim=(-numpy.pi, +numpy.pi),
    xlabel="Repetition", ylabel="Phase (rad)")
plots[0].legend()

figure.savefig(destination)

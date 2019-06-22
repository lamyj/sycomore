import math
import os
import struct
import unittest

import sycomore
from sycomore.units import *

class TestOffResonance(unittest.TestCase):
    def test_off_resonance(self):
        species = sycomore.Species(0*Hz, 0*Hz, 0*um*um/ms)
        m0 = sycomore.Magnetization(0, 0, 1)

        pulse = sycomore.Pulse(90*deg, math.pi*rad)
        pulse_duration = 1*ms
        pulse_support_size = 101
        zero_crossings = 2

        # NOTE: in the absence of relaxation and diffusion, the TR is meaningless
        TR = 500*ms;
        slice_thickness = 1*mm;

        t0 = pulse_duration/(2*zero_crossings)
        sinc_pulse = sycomore.HardPulseApproximation(
            pulse,
            sycomore.linspace(pulse_duration, pulse_support_size),
            sycomore.sinc_envelope(t0), 1/t0, slice_thickness, "rf")

        refocalization = sycomore.TimeInterval(
            (TR-pulse_duration)/2., -sinc_pulse.get_gradient_moment()/2)

        model = sycomore.como.Model(
            species, m0, [
                ["rf", sinc_pulse.get_time_interval()],
                ["refocalization", refocalization]])

        model.apply_pulse(sinc_pulse)
        model.apply_time_interval("refocalization")

        frequencies = sycomore.linspace(60.*rad/ms, 201)
        magnetization = [
            model.isochromat(set(), sycomore.Point(), f) for f in frequencies]

        root = os.environ["SYCOMORE_TEST_DATA"]
        with open(os.path.join(root, "baseline", "off_resonance.dat"), "rb") as fd:
            contents = fd.read()
            baseline = struct.unpack((int(len(contents)/8))*"d", contents)

        self.assertEqual(len(baseline), 2*len(magnetization))
        for i in range(len(magnetization)):
            self.assertAlmostEqual(
                sycomore.transversal(magnetization[i]), baseline[2*i])
            self.assertAlmostEqual(magnetization[i][2], baseline[2*i+1])

if __name__ == "__main__":
    unittest.main()

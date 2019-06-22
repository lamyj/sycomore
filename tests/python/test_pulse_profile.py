import math
import os
import struct
import unittest

import sycomore
from sycomore.units import *

class TestPulseProfile(unittest.TestCase):
    def test_pulse_profile(self):
        species = sycomore.Species(0*Hz, 0*Hz, 0*um*um/ms)
        m0 = sycomore.Magnetization(0, 0, 1)

        pulse = sycomore.Pulse(90*deg, math.pi*rad)
        pulse_duration = 1*ms
        pulse_support_size = 101
        zero_crossings = 2

        # NOTE: in the absence of relaxation and diffusion, the TR is meaningless
        TR = 500*ms;
        slice_thickness = 1*mm;

        sampling_support_size = 501

        t0 = pulse_duration/(2*zero_crossings)
        sinc_pulse = sycomore.HardPulseApproximation(
            pulse,
            sycomore.linspace(pulse_duration, pulse_support_size),
            sycomore.sinc_envelope(t0), 1/t0, slice_thickness, "rf")

        refocalization = sycomore.TimeInterval(
            (TR-pulse_duration)/2., -sinc_pulse.get_gradient_moment()/2)

        sampling_locations = sycomore.linspace(
            sycomore.Point(0*m, 0*m, 2*slice_thickness), sampling_support_size)

        model = sycomore.como.Model(
            species, m0, [
                ["rf", sinc_pulse.get_time_interval()],
                ["refocalization", refocalization]])

        model.apply_pulse(sinc_pulse)

        before_refocalization = [
            model.isochromat(set(), p) for p in sampling_locations]

        model.apply_time_interval("refocalization")

        after_refocalization = [
            model.isochromat(set(), p) for p in sampling_locations]

        root = os.environ["SYCOMORE_TEST_DATA"]
        with open(os.path.join(root, "baseline", "pulse_profile.dat"), "rb") as fd:
            contents = fd.read()
            baseline = struct.unpack((int(len(contents)/8))*"d", contents)

        self.assertEqual(len(baseline), 2*3*len(sampling_locations))
        for i in range(len(sampling_locations)):
            m_test = before_refocalization[i]
            m_baseline = baseline[3*i:3*(i+1)]

            self.assertAlmostEqual(m_test[0], m_baseline[0])
            self.assertAlmostEqual(m_test[1], m_baseline[1])
            self.assertAlmostEqual(m_test[2], m_baseline[2])
        for i in range(len(sampling_locations)):
            m_test = after_refocalization[i]
            m_baseline = baseline[
                3*(i+len(sampling_locations)):3*(i+len(sampling_locations)+1)]

            self.assertAlmostEqual(m_test[0], m_baseline[0])
            self.assertAlmostEqual(m_test[1], m_baseline[1])
            self.assertAlmostEqual(m_test[2], m_baseline[2])

if __name__ == "__main__":
    unittest.main()

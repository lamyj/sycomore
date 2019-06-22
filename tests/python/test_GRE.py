import math
import os
import struct
import unittest

import sycomore
from sycomore.units import *

class TestGRE(unittest.TestCase):
    def setUp(self):
        self.species = sycomore.Species(1000*ms, 100*ms, 0.89*um*um/ms)
        self.m0 = sycomore.Magnetization(0, 0, 1)

        self.flip_angle = 40*deg
        self.pulse_duration = 1*ms
        self.pulse_support_size = 101
        self.zero_crossings = 2

        self.TR = 500*ms
        self.slice_thickness = 1*mm

        self.TR_count = 10

    def test_ideal(self):
        model = sycomore.como.Model(
            self.species, self.m0,
            [["half_echo", sycomore.TimeInterval(self.TR/2.)]])

        magnetization = []
        for i in range(self.TR_count):
            model.apply_pulse(
                sycomore.Pulse(self.flip_angle, (math.pi/3+(i%2)*math.pi)*rad))
            model.apply_time_interval("half_echo")
            magnetization.append(model.isochromat())
            model.apply_time_interval("half_echo")

        root = os.environ["SYCOMORE_TEST_DATA"]
        with open(os.path.join(root, "baseline", "GRE_ideal.dat"), "rb") as fd:
            contents = fd.read()
            baseline = struct.unpack((int(len(contents)/8))*"d", contents)

        self.assertEqual(len(baseline), 3*self.TR_count)
        for i in range(self.TR_count):
            m_test = magnetization[i]
            m_baseline = baseline[3*i:3*(i+1)]

            self.assertAlmostEqual(m_test[0], m_baseline[0])
            self.assertAlmostEqual(m_test[1], m_baseline[1])
            self.assertAlmostEqual(m_test[2], m_baseline[2])

    def test_real(self):
        t0 = self.pulse_duration/(2*self.zero_crossings)
        sinc_pulse = sycomore.HardPulseApproximation(
            sycomore.Pulse(self.flip_angle, 0*rad),
            sycomore.linspace(self.pulse_duration, self.pulse_support_size),
            sycomore.sinc_envelope(t0), 1/t0, self.slice_thickness, "rf")

        half_echo = sycomore.TimeInterval(
            (self.TR-self.pulse_duration)/2.,
            -sinc_pulse.get_gradient_moment()/2)

        model = sycomore.como.Model(
            self.species, self.m0, [
                ["rf", sinc_pulse.get_time_interval()],
                ["half_echo", half_echo]])

        magnetization = []
        for i in range(self.TR_count):
            sinc_pulse.set_phase((math.pi/3+(i%2)*math.pi)*rad)
            model.apply_pulse(sinc_pulse)
            model.apply_time_interval("half_echo")
            magnetization.append(model.isochromat())
            model.apply_time_interval("half_echo")

        root = os.environ["SYCOMORE_TEST_DATA"]
        with open(os.path.join(root, "baseline", "GRE_real.dat"), "rb") as fd:
            contents = fd.read()
            baseline = struct.unpack((int(len(contents)/8))*"d", contents)

        self.assertEqual(len(baseline), 3*self.TR_count)
        for i in range(self.TR_count):
            m_test = magnetization[i]
            m_baseline = baseline[3*i:3*(i+1)]

            self.assertAlmostEqual(m_test[0], m_baseline[0])
            self.assertAlmostEqual(m_test[1], m_baseline[1])
            self.assertAlmostEqual(m_test[2], m_baseline[2])

if __name__ == "__main__":
    unittest.main()

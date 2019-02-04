import math
import unittest

import sycomore
from sycomore.units import *

class TestPulse(unittest.TestCase):
    def test_constructor_real(self):
        pulse = sycomore.Pulse(1,2)
        self.assertEqual(pulse.angle, 1)
        self.assertEqual(pulse.phase, 2)

    def test_constructor_quantities(self):
        pulse = sycomore.Pulse(1*deg,2*deg)
        self.assertEqual(pulse.angle, 1/180*math.pi)
        self.assertEqual(pulse.phase, 2/180*math.pi)
        with self.assertRaises(Exception):
            sycomore.Pulse(m, rad)
        with self.assertRaises(Exception):
            sycomore.Pulse(rad, m)

if __name__ == "__main__":
    unittest.main()

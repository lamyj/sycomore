import math
import unittest

import sycomore
from sycomore.units import *

class TestPulse(unittest.TestCase):
    def test_constructor_quantities(self):
        pulse = sycomore.Pulse(1*deg,2*deg)
        self.assertEqual(pulse.angle, 1*deg)
        self.assertEqual(pulse.phase, 2*deg)
        with self.assertRaises(Exception):
            sycomore.Pulse(m, rad)
        with self.assertRaises(Exception):
            sycomore.Pulse(rad, m)

if __name__ == "__main__":
    unittest.main()

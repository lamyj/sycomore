import math
import unittest

import sycomore
from sycomore.units import *

class TestPulse(unittest.TestCase):
    def test_constructor_quantities(self):
        pulse_1 = sycomore.Pulse(1*deg,2*deg)
        self.assertEqual(pulse_1.angle, 1*deg)
        self.assertEqual(pulse_1.phase, 2*deg)
        
        pulse_2 = sycomore.Pulse(2*deg)
        self.assertEqual(pulse_2.angle, 2*deg)
        self.assertEqual(pulse_2.phase, 0*deg)
        
        with self.assertRaises(Exception):
            sycomore.Pulse(m, rad)
        with self.assertRaises(Exception):
            sycomore.Pulse(rad, m)

if __name__ == "__main__":
    unittest.main()

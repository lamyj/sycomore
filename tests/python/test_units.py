import unittest

import sycomore
from sycomore.units import *

class TestUnits(unittest.TestCase):
    def test_basic_unit(self):
        length = 100*cm
        self.assertEqual(length.magnitude, 1)
        self.assertEqual(length.dimensions, sycomore.Length)

    def test_derived_unit(self):
        length = 1*kN
        self.assertEqual(length.magnitude, 1000)
        self.assertEqual(length.dimensions, sycomore.Force)

if __name__ == "__main__":
    unittest.main()

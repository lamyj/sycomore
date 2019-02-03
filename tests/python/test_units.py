import unittest

import sycomore

class TestUnits(unittest.TestCase):
    def test_basic_unit(self):
        length = 100*sycomore.units.cm
        self.assertEqual(length.magnitude, 1)
        self.assertEqual(length.dimensions, sycomore.Length)

    def test_derived_unit(self):
        length = 1*sycomore.units.kN
        self.assertEqual(length.magnitude, 1000)
        self.assertEqual(length.dimensions, sycomore.Force)

if __name__ == "__main__":
    unittest.main()

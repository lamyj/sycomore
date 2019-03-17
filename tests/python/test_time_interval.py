import unittest

import sycomore
from sycomore.units import *

class TestTimeInterval(unittest.TestCase):
    def test_quantity_scalar_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, 2./dm)
        self.assertEqual(interval.duration, 1e-3*s)
        self.assertSequenceEqual(
            interval.gradient_moment, [20.*rad/m,20.*rad/m,20.*rad/m])

    def test_quantity_vector_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, [2*rad/dm, 4*rad/m, 8*rad/dam])
        self.assertEqual(interval.duration, 1e-3*s)
        self.assertSequenceEqual(
            interval.gradient_moment, [20*rad/m, 4*rad/m, 0.8*rad/m])

if __name__ == "__main__":
    unittest.main()

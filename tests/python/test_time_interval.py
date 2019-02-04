import unittest

import sycomore
from sycomore.units import *

class TestTimeInterval(unittest.TestCase):
    def test_scalar_scalar_constructor(self):
        interval = sycomore.TimeInterval(1., 2.)
        self.assertEqual(interval.duration, 1.)
        self.assertSequenceEqual(interval.gradient_moment, [2.,2.,2.])

    def test_quantity_scalar_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, 2./dm)
        self.assertEqual(interval.duration, 1e-3)
        self.assertSequenceEqual(interval.gradient_moment, [20.,20.,20.])

    def test_scalar_array_constructor(self):
        interval = sycomore.TimeInterval(1., [2., 3., 4.])
        self.assertEqual(interval.duration, 1.)
        self.assertSequenceEqual(interval.gradient_moment, [2., 3., 4.])

    def test_quantity_vector_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, [2/dm, 4/m, 8/dam])
        self.assertEqual(interval.duration, 1e-3)
        self.assertSequenceEqual(interval.gradient_moment, [20, 4, 0.8])

if __name__ == "__main__":
    unittest.main()

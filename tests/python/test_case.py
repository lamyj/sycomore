import unittest
import numpy

import sycomore

def assertQuantityAlmostEqual(left, right, msg=None):
    if not numpy.allclose(left.magnitude, right.magnitude):
        raise AssertionError(
            msg if msg is not None
            else f"Unequal magnitudes: {left.magnitude} vs. {right.magnitude}")
    if left.dimensions != right.dimensions:
        raise AssertionError(
            msg if msg is not None
            else f"Unequal dimensions: {left.dimensions} vs. {right.dimensions}")

class TestCase(unittest.TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(sycomore.Matrix3x3Q, assertQuantityAlmostEqual)
        self.addTypeEqualityFunc(sycomore.ArrayQ, assertQuantityAlmostEqual)

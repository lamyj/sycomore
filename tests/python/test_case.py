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
        self.addTypeEqualityFunc(sycomore.Vector3Q, assertQuantityAlmostEqual)
        self.addTypeEqualityFunc(sycomore.Matrix3x3Q, assertQuantityAlmostEqual)
        self.addTypeEqualityFunc(sycomore.ArrayQ, assertQuantityAlmostEqual)
    
    def assertAllClose(self, left, right, **kwargs):
        left_dimensions = (
            left.dimensions if hasattr(left, "dimensions")
            else sycomore.Dimensions())
        right_dimensions = (
            right.dimensions if hasattr(right, "dimensions")
            else sycomore.Dimensions())
        if left_dimensions != right_dimensions:
            raise AssertionError(
                f"Unequal dimensions: "
                f"{left_dimensions} vs. {right_dimensions}")
        
        left_magnitude = left.magnitude if hasattr(left, "magnitude") else left
        right_magnitude = right.magnitude if hasattr(right, "magnitude") else right
        return numpy.testing.assert_allclose(
            left_magnitude, right_magnitude, **kwargs)

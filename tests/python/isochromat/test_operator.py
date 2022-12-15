import numpy
import unittest

import sycomore

class TestOperator(unittest.TestCase):
    def setUp(self):
        self.left = numpy.array([
            [
                [1, 2, 3, 4],
                [5, 6, 7, 8],
                [9, 10, 11, 12],
                [13, 14, 15, 16]
            ], [
                [13, 9, 5, 1],
                [14, 10, 6, 2],
                [15, 11, 7, 3],
                [16, 12, 8, 4]
            ]])
        
        self.right = [
            [
                [16, 15, 14, 13],
                [12, 11, 10, 9],
                [8, 7, 6, 5],
                [4, 3, 2, 1]
            ], [
                [4, 3, 2, 1],
                [8, 7, 6, 5],
                [12, 11, 10, 9],
                [16, 15, 14, 13]
            ]]
        
    def test_empty_constructor(self):
        operator = sycomore.isochromat.Operator()
        numpy.testing.assert_almost_equal(
            operator.array, numpy.identity(4).reshape(1, 4, 4))
    
    def test_data_constructor(self):
        operator = sycomore.isochromat.Operator(self.left)
        numpy.testing.assert_almost_equal(operator.array, self.left)
    
    def test_inplace_product(self):
        left = sycomore.isochromat.Operator(self.left)
        right = sycomore.isochromat.Operator(self.right)
        
        right.preMultiply(left)
    
        numpy.testing.assert_almost_equal(right.array, self.left @ self.right)
    
    def test_pre_multiply(self):
        left = sycomore.isochromat.Operator(self.left)
        right = sycomore.isochromat.Operator(self.right)
        
        left *= right
    
        numpy.testing.assert_almost_equal(left.array, self.left @ self.right)
    
    def test_multiply(self):
        left = sycomore.isochromat.Operator(self.left)
        right = sycomore.isochromat.Operator(self.right)
        
        combined = left * right
    
        numpy.testing.assert_almost_equal(
            combined.array, self.left @ self.right)
if __name__ == "__main__":
    unittest.main()

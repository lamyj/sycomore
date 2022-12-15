import numpy
import unittest

import sycomore
from sycomore.units import *

class TestArray(unittest.TestCase):
    def test_empty_constructor(self):
        array = sycomore.Array[float]()
        self.assertEqual(array.size(), 0)
        self.assertTrue(array.empty())
        self.assertFalse(array)

    def test_buffer_constructor(self):
        array = sycomore.Array[float](numpy.array([1.,2.,3.]))
        self.assertEqual(array.size(), 3)
        self.assertFalse(array.empty())
        self.assertTrue(array)

    def test_sequence_constructor(self):
        array = sycomore.Array[float]([1,2,3])
        self.assertEqual(array.size(), 3)
        self.assertFalse(array.empty())
        self.assertTrue(array)

    def test_args_constructor(self):
        array = sycomore.Array[float](1,2,3)
        self.assertEqual(array.size(), 3)
        self.assertFalse(array.empty())
        self.assertTrue(array)

    def test_item(self):
        array = sycomore.Array[float]([1,2,3])
        array[1] = 42
        self.assertEqual(array[1], 42)

    def test_iter(self):
        array = sycomore.Array[float]([1,2,3])
        self.assertSequenceEqual(array, [1,2,3])

    def test_buffer(self):
        array = sycomore.Array[float]([1,2,3])
        self.assertTrue((numpy.array(array, copy=False) == [1,2,3]).all())
    
    def test_operators(self):
        array = sycomore.Array[float]([1,2,3])
        self.assertTrue(array+4. == sycomore.Array[float]([5,6,7]))
        self.assertTrue(array-4. == sycomore.Array[float]([-3,-2,-1]))
        self.assertTrue(array*4. == sycomore.Array[float]([4,8,12]))
        self.assertTrue(array/4. == sycomore.Array[float]([0.25,0.5,0.75]))
        
        self.assertTrue(4+array == sycomore.Array[float]([5,6,7]))
        self.assertTrue(4-array == sycomore.Array[float]([3,2,1]))
        self.assertTrue(4*array == sycomore.Array[float]([4,8,12]))
        self.assertTrue(3/array == sycomore.Array[float]([3,1.5,1]))
if __name__ == "__main__":
    unittest.main()

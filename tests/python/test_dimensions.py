import unittest

import sycomore

class TestDimensions(unittest.TestCase):
    def test_comparison(self):
        self.assertTrue(sycomore.Length == sycomore.Length)
        self.assertTrue(sycomore.Length != sycomore.Time)

    def test_multiplication_in_place(self):
        d1 = sycomore.Dimensions(1,0,0,0,0,0,0)
        d2 = sycomore.Dimensions(0,1,0,0,0,0,0)
        r = sycomore.Dimensions(1,1,0,0,0,0,0)
        d1 *= d2
        self.assertEqual(d1, r)

    def test_division_in_place(self):
        d1 = sycomore.Dimensions(1,0,0,0,0,0,0)
        d2 = sycomore.Dimensions(0,1,0,0,0,0,0)
        r = sycomore.Dimensions(1,-1,0,0,0,0,0)
        d1 /= d2
        self.assertEqual(d1, r)

    def test_multiplication(self):
        d1 = sycomore.Dimensions(1,0,0,0,0,0,0)
        d2 = sycomore.Dimensions(0,1,0,0,0,0,0)
        r = sycomore.Dimensions(1,1,0,0,0,0,0)
        self.assertEqual(d1*d2, r)

    def test_division(self):
        d1 = sycomore.Dimensions(1,0,0,0,0,0,0)
        d2 = sycomore.Dimensions(0,1,0,0,0,0,0)
        r = sycomore.Dimensions(1,-1,0,0,0,0,0)
        self.assertEqual(d1/d2, r)

    def test_pow(self):
        d = sycomore.Dimensions(3,0,-2,0,0,0,0)
        r = sycomore.Dimensions(6,0,-4,0,0,0,0)
        self.assertEqual(d**2, r)

if __name__ == "__main__":
    unittest.main()

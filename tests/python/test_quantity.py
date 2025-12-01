import logging
import math
import sys
import pickle
import unittest

import numpy

import sycomore

class TestQuantity(unittest.TestCase):
    def test_comparison(self):
        q1 = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        q3 = sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0))
        scalar = sycomore.Quantity(2, sycomore.Dimensions(0,0,0,0,0,0,0))
        
        self.assertTrue(q1 == q1)
        self.assertFalse(q1 == q2)
        self.assertFalse(q1 == q3)
        self.assertTrue(scalar == 2)
        self.assertTrue(2 == scalar)
        self.assertFalse(scalar == 3)
        self.assertFalse(3 == scalar)

        self.assertFalse(q1 != q1)
        self.assertTrue(q1 != q2)
        self.assertTrue(q1 != q3)
        self.assertFalse(scalar != 2)
        self.assertFalse(2 != scalar)
        self.assertTrue(scalar != 3)
        self.assertTrue(3 !=scalar)

        self.assertTrue(q1 < q2)
        self.assertTrue(scalar < 3)
        self.assertTrue(1 < scalar)
        
        self.assertFalse(q2 <= q1)
        self.assertTrue(scalar <= 2)
        self.assertTrue(2 <= scalar)
        
        self.assertFalse(q1 > q2)
        self.assertTrue(scalar > 1)
        self.assertTrue(3 > scalar)
        
        self.assertFalse(q1 >= q2)
        self.assertTrue(scalar >= 2)
        self.assertTrue(2 >= scalar)
        
        with self.assertRaises(Exception):
            q1 < q3
        with self.assertRaises(Exception):
            q1 <= q3
        with self.assertRaises(Exception):
            q1 > q3
        with self.assertRaises(Exception):
            q1 >= q3
        
        with self.assertRaises(Exception):
            q1 < scalar
        with self.assertRaises(Exception):
            scalar < q1
        with self.assertRaises(Exception):
            q1 <= scalar
        with self.assertRaises(Exception):
            scalar <= q1
        with self.assertRaises(Exception):
            q1 > scalar
        with self.assertRaises(Exception):
            scalar > q1
        with self.assertRaises(Exception):
            q1 >= scalar
        with self.assertRaises(Exception):
            scalar >= q1

    def test_addition_in_place(self):
        q1 = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(5, sycomore.Dimensions(1,0,0,0,0,0,0))
        q1 += q2
        self.assertEqual(q1, r1)
        
        scalar = sycomore.Quantity(2, sycomore.Dimensions(0,0,0,0,0,0,0))
        r2 = sycomore.Quantity(5, sycomore.Dimensions(0,0,0,0,0,0,0))
        scalar += 3
        self.assertEqual(scalar, r2)

        q3 = sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0))
        with self.assertRaises(Exception):
            q1 += q3
        with self.assertRaises(Exception):
            q1 += 3

    def test_subtraction_in_place(self):
        q1 = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(-1, sycomore.Dimensions(1,0,0,0,0,0,0))
        q1 -= q2
        self.assertEqual(q1, r1)
        
        scalar = sycomore.Quantity(2, sycomore.Dimensions(0,0,0,0,0,0,0))
        r2 = sycomore.Quantity(-1, sycomore.Dimensions(0,0,0,0,0,0,0))
        scalar -= 3
        self.assertEqual(scalar, r2)

        q3 = sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0))
        with self.assertRaises(Exception):
            q1 -= q3
        with self.assertRaises(Exception):
            q1 -= 3

    def test_multiplication_in_place(self):
        q1 = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(0,1,0,0,0,0,0))
        r = sycomore.Quantity(6, sycomore.Dimensions(1,1,0,0,0,0,0))
        q1 *= q2
        self.assertEqual(q1, r)

    def test_scalar_multiplication_in_place(self):
        q = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(6, sycomore.Dimensions(1,0,0,0,0,0,0))
        q *= 3
        self.assertEqual(q, r)

    def test_division_in_place(self):
        q1 = sycomore.Quantity(1, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(4, sycomore.Dimensions(0,1,0,0,0,0,0))
        r = sycomore.Quantity(0.25, sycomore.Dimensions(1,-1,0,0,0,0,0))
        q1 /= q2
        self.assertEqual(q1, r)

    def test_scalar_division_in_place(self):
        q = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(0.5, sycomore.Dimensions(1,0,0,0,0,0,0))
        q /= 4
        self.assertEqual(q, r)
    
    def test_floordiv_in_place(self):
        q1 = sycomore.Quantity(7, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(0,1,0,0,0,0,0))
        r = sycomore.Quantity(2, sycomore.Dimensions(1,-1,0,0,0,0,0))
        q1 //= q2
        self.assertEqual(q1, r)

    def test_scalar_floordiv_in_place(self):
        q = sycomore.Quantity(7, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(7//3, sycomore.Dimensions(1,0,0,0,0,0,0))
        q //= 3
        self.assertEqual(q, r)

    def test_modulo_in_place(self):
        q1 = sycomore.Quantity(7, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(1, sycomore.Dimensions(1,0,0,0,0,0,0))
        q1 %= q2
        self.assertEqual(q1, r)

        q3 = sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0))
        with self.assertRaises(Exception):
            q1 %= q3

    def test_scalar_modulo_in_place(self):
        q = sycomore.Quantity(7, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(1, sycomore.Dimensions(1,0,0,0,0,0,0))
        q %= 3
        self.assertEqual(q, r)
    
    def test_convert_to(self):
        q1 = sycomore.Quantity(70, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(10, sycomore.Dimensions(1,0,0,0,0,0,0))
        q3 = sycomore.Quantity(10, sycomore.Dimensions(0,1,0,0,0,0,0))
        r = 7
        self.assertEqual(q1.convert_to(q2), r)
        with self.assertRaises(Exception):
            q1.convert_to(q3)
    
    def test_float(self):
        scalar = sycomore.Quantity(3, sycomore.Dimensions(0,0,0,0,0,0,0))
        q = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(float(scalar), 3)
        with self.assertRaises(Exception):
            float(q)

    def test_unary_plus(self):
        q = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(+q, q)

    def test_unary_minus(self):
        q = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(-2, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(-q, r)

    def test_addition(self):
        q1 = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(5, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(q1+q2, r1)
        
        scalar = sycomore.Quantity(2, sycomore.Dimensions(0,0,0,0,0,0,0))
        r2 = sycomore.Quantity(5, sycomore.Dimensions(0,0,0,0,0,0,0))
        self.assertEqual(scalar+3, r2)
        self.assertEqual(3+scalar, r2)
        
        q3 = sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0))
        with self.assertRaises(Exception):
            q1 + q3
        with self.assertRaises(Exception):
            q1 + 3
        with self.assertRaises(Exception):
            3 + q1

    def test_subtraction(self):
        q1 = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(-1, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(q1-q2, r1)
        
        scalar = sycomore.Quantity(2, sycomore.Dimensions(0,0,0,0,0,0,0))
        r2 = sycomore.Quantity(-1, sycomore.Dimensions(0,0,0,0,0,0,0))
        self.assertEqual(scalar-3, r2)
        self.assertEqual(1-scalar, r2)

        q3 = sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0))
        with self.assertRaises(Exception):
            q1 - q3
        with self.assertRaises(Exception):
            q1 - 3
        with self.assertRaises(Exception):
            1 - q1

    def test_multiplication(self):
        q1 = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(0,1,0,0,0,0,0))
        r = sycomore.Quantity(6, sycomore.Dimensions(1,1,0,0,0,0,0))
        self.assertEqual(q1*q2, r)

    def test_scalar_multiplication(self):
        q = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(6, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(q*3, r)
        self.assertEqual(3*q, r)

    def test_division(self):
        q1 = sycomore.Quantity(1, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(4, sycomore.Dimensions(0,1,0,0,0,0,0))
        r = sycomore.Quantity(0.25, sycomore.Dimensions(1,-1,0,0,0,0,0))
        self.assertEqual(q1/q2, r)

    def test_scalar_division(self):
        q = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(0.5, sycomore.Dimensions(1,0,0,0,0,0,0))
        r2 = sycomore.Quantity(1.5, sycomore.Dimensions(-1,0,0,0,0,0,0))
        self.assertEqual(q/4, r1)
        self.assertEqual(3/q, r2)

    def test_floordiv(self):
        q1 = sycomore.Quantity(7, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(0,1,0,0,0,0,0))
        r = sycomore.Quantity(7//3, sycomore.Dimensions(1,-1,0,0,0,0,0))
        self.assertEqual(q1//q2, r)
    
    def test_floordiv_scalar(self):
        q = sycomore.Quantity(7, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(7//3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r2 = sycomore.Quantity(15//7, sycomore.Dimensions(-1,0,0,0,0,0,0))
        self.assertEqual(q//3, r1)
        self.assertEqual(15//q, r2)

    def test_modulo(self):
        q1 = sycomore.Quantity(7, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(1, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(q1%q2, r)

    def test_scalar_modulo(self):
        q = sycomore.Quantity(7, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(1, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(q%3, r)
    
    def test_divmod(self):
        q1 = sycomore.Quantity(7, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = (
            sycomore.Quantity(2, sycomore.Dimensions(0,0,0,0,0,0,0)),
            sycomore.Quantity(1, sycomore.Dimensions(1,0,0,0,0,0,0)))
        self.assertEqual(divmod(q1, q2), r)
    
    def test_divmod_scalar(self):
        q = sycomore.Quantity(7, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = (
            sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0)),
            sycomore.Quantity(1, sycomore.Dimensions(1,0,0,0,0,0,0)))
        self.assertEqual(divmod(q, 3), r)
    
    def test_abs(self):
        q1 = sycomore.Quantity(-9, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(9, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(abs(q1), r1)
        
        q2 = sycomore.Quantity(9, sycomore.Dimensions(-1,0,0,0,0,0,0))
        r2 = sycomore.Quantity(9, sycomore.Dimensions(-1,0,0,0,0,0,0))
        self.assertEqual(abs(q2), r2)

    def test_pow(self):
        q = sycomore.Quantity(9, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(3, sycomore.Dimensions(0.5,0,0,0,0,0,0))
        self.assertEqual(q**0.5, r)

    def test_round(self):
        q1 = sycomore.Quantity(9.2, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(9, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(round(q1), r1)
        
        q2 = sycomore.Quantity(-9.7, sycomore.Dimensions(1,0,0,0,0,0,0))
        r2 = sycomore.Quantity(-10, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(round(q2), r2)
        
        q3 = sycomore.Quantity(9, sycomore.Dimensions(-1.5,0,0,0,0,0,0))
        r3 = sycomore.Quantity(9, sycomore.Dimensions(-1.5,0,0,0,0,0,0))
        self.assertEqual(round(q3), r3)

    def test_trunc(self):
        q1 = sycomore.Quantity(9.2, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(9, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(math.trunc(q1), r1)
        
        q2 = sycomore.Quantity(-9.7, sycomore.Dimensions(1,0,0,0,0,0,0))
        r2 = sycomore.Quantity(-9, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(math.trunc(q2), r2)
        
        q3 = sycomore.Quantity(9, sycomore.Dimensions(-1.5,0,0,0,0,0,0))
        r3 = sycomore.Quantity(9, sycomore.Dimensions(-1.5,0,0,0,0,0,0))
        self.assertEqual(math.trunc(q3), r3)
    
    def test_floor(self):
        q1 = sycomore.Quantity(9.2, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(9, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(math.floor(q1), r1)
        
        q2 = sycomore.Quantity(-9.7, sycomore.Dimensions(1,0,0,0,0,0,0))
        r2 = sycomore.Quantity(-10, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(math.floor(q2), r2)
        
        q3 = sycomore.Quantity(9, sycomore.Dimensions(-1.5,0,0,0,0,0,0))
        r3 = sycomore.Quantity(9, sycomore.Dimensions(-1.5,0,0,0,0,0,0))
        self.assertEqual(math.floor(q3), r3)
    
    def test_ceil(self):
        q1 = sycomore.Quantity(9.2, sycomore.Dimensions(1,0,0,0,0,0,0))
        r1 = sycomore.Quantity(10, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(math.ceil(q1), r1)
        
        q2 = sycomore.Quantity(-9.7, sycomore.Dimensions(1,0,0,0,0,0,0))
        r2 = sycomore.Quantity(-9, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(math.ceil(q2), r2)
        
        q3 = sycomore.Quantity(9, sycomore.Dimensions(-1.5,0,0,0,0,0,0))
        r3 = sycomore.Quantity(9, sycomore.Dimensions(-1.5,0,0,0,0,0,0))
        self.assertEqual(math.ceil(q3), r3)

    def test_pickle(self):
        q = sycomore.Quantity(0.5, sycomore.Dimensions(7,6,5,4,3,2,1))
        self.assertEqual(pickle.loads(pickle.dumps(q)), q)

    def test_hash(self):
        quantities = set()

        quantities.add(sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0)))
        quantities.add(sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0)))
        quantities.add(sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0)))
        self.assertEqual(len(quantities), 3)
        
        quantities.add(sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0)))
        self.assertEqual(len(quantities), 3)

if __name__ == "__main__":
    unittest.main()

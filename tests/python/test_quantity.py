import logging
import math
import sys
if sys.version_info[0] == 2:
    import cPickle
else:
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
        if sys.version_info[0] == 2:
            # WARNING: when running with Python2, only cPickle with version >= 2
            # works. Refer to the last paragraph of 
            # https://pybind11.readthedocs.io/en/stable/advanced/classes.html?highlight=pickle#pickling-support
            self.assertEqual(cPickle.loads(cPickle.dumps(q, -1)), q)
        else:
            self.assertEqual(pickle.loads(pickle.dumps(q)), q)

    def test_hash(self):
        quantities = set()

        quantities.add(sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0)))
        quantities.add(sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0)))
        quantities.add(sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0)))
        self.assertEqual(len(quantities), 3)
        
        quantities.add(sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0)))
        self.assertEqual(len(quantities), 3)
    
    def test_ufuncs(self):
        ufuncs = [
            x for x in dir(numpy) if isinstance(getattr(numpy, x), numpy.ufunc)]
        non_object_ufuncs = [
            x for x in ufuncs if not any(
                t.endswith("O->O") for t in getattr(numpy, x).types)
        ]
        logging.warning(
            "The following ufuncs do not operate on objects: {}".format(
                ", ".join(non_object_ufuncs)))
        ufuncs = [x for x in ufuncs if x not in non_object_ufuncs]
        
        not_applicable_ufuncs = [
            "conj", "conjugate", # Quantities are real-valued
            "deg2rad", "degrees", "rad2deg", "radians", 
            "bitwise_and", "bitwise_or", "bitwise_xor", "bitwise_not", "invert", 
            "left_shift", "right_shift", "logical_and", "logical_or", 
            "logical_xor", "logical_not", "isnat", 
            "gcd", "lcm", "matmul"
        ]
        ufuncs = [x for x in ufuncs if x not in not_applicable_ufuncs]
        
        from sycomore.units import m, deg
        Scalar = sycomore.Quantity(1, sycomore.Dimensions())
        tests = [
            ["add", [1*m, 2*m], 3*m],
            ["add", [1*Scalar, 2], 3*Scalar],
            ["add", [1, 2*Scalar], 3*Scalar],
            
            ["subtract", [1*m, 2*m], -1*m],
            ["subtract", [1*Scalar, 2], -1*Scalar],
            ["subtract", [1, 2*Scalar], -1*Scalar],
            
            ["multiply", [2*m, 3*m], 6*m**2],
            ["multiply", [2*m, 3], 6*m],
            ["multiply", [2, 3*m], 6*m],
            
            ["divide", [3*m, 2*m], 1.5*Scalar],
            ["divide", [3*m, 2], 1.5*m],
            ["divide", [3, 2*m], 1.5/m],
            # Same tests for true_divide
            
            ["floor_divide", [3*m, 2*m], 1*Scalar],
            ["floor_divide", [3*m, 2], 1*m],
            ["floor_divide", [3, 2*m], 1/m],
            
            ["negative", [-1*m], 1*m],
            ["positive", [-1*m], -1*m],
            ["power", [3*m, 2], 9*m**2],
            
            ["remainder", [3*m, 2], 1*m],
            ["remainder", [3*m, 2*m], 1*m],
            # Same tests for mod and fmod
            
            ["absolute", [-1*m], 1*m],
            # Same tests for fabs and abs
            
            ["rint", [-1.7*m], -2*m],
            
            ["sign", [-1.7*Scalar], -1],
            
            ["exp", [0*Scalar], 1*Scalar],
            ["exp2", [1*Scalar], 2*Scalar],
            ["log", [1*Scalar], 0*Scalar],
            ["log2", [8*Scalar], 3*Scalar],
            ["log10", [1000*Scalar], 3*Scalar],
            ["expm1", [0*Scalar], 0*Scalar],
            ["log1p", [0*Scalar], 0*Scalar],
            
            ["sqrt", [4*m**2], 2*m],
            ["square", [2*m], 4*m**2],
            ["cbrt", [8*m**3], 2*m],
            ["reciprocal", [2*m], 0.5/m],
            
            ["sin", [0*deg], 0*Scalar],
            ["cos", [0*deg], 1*Scalar],
            ["tan", [0*deg], 0*Scalar],
            ["arcsin", [0*Scalar], 0*deg],
            ["arccos", [1*Scalar], 0*deg],
            ["arctan", [0*Scalar], 0*deg],
            ["arctan2", [0*Scalar, 1*Scalar], 0*deg],
            ["hypot", [3*Scalar, 4*Scalar], 5*Scalar],
            
            ["sinh", [0*deg], 0*Scalar],
            ["cosh", [0*deg], 1*Scalar],
            ["tanh", [0*deg], 0*Scalar],
            ["arcsinh", [0*Scalar], 0*deg],
            ["arccosh", [1*Scalar], 0*deg],
            ["arctanh", [0*Scalar], 0*deg],
            
            ["greater", [2*m, 1*m], True],
            ["greater", [2*Scalar, 1], True],
            ["greater", [2, 1*Scalar], True],
            
            ["greater_equal", [2*m, 2*m], True],
            ["greater_equal", [2*Scalar, 2], True],
            ["greater_equal", [2, 2*Scalar], True],
            
            ["less", [1*m, 2*m], True],
            ["less", [1*Scalar, 2], True],
            ["less", [1, 2*Scalar], True],
            
            ["less_equal", [2*m, 2*m], True],
            ["less_equal", [2*Scalar, 2], True],
            ["less_equal", [2, 2*Scalar], True],
            
            ["not_equal", [2*m, 1*m], True],
            ["not_equal", [2*Scalar, 1], True],
            ["not_equal", [2, 1*Scalar], True],
            
            ["equal", [2*m, 2*m], True],
            ["equal", [2*Scalar, 2], True],
            ["equal", [2, 2*Scalar], True],
            
            ["maximum", [1*m, 2*m], 2*m],
            ["maximum", [1*Scalar, 2], 2*Scalar],
            ["maximum", [1, 2*Scalar], 2*Scalar],
            
            ["minimum", [1*m, 2*m], 1*m],
            ["minimum", [1*Scalar, 2], 1*Scalar],
            ["minimum", [1, 2*Scalar], 1*Scalar],
            
            ["ceil", [-2.1*m], -2*m],
            ["floor", [-2.1*m], -3*m],
            ["trunc", [-2.1*m], -2*m],
        ]
        equivalences = [
            ["true_divide", "divide"], ["mod", "remainder"], 
            ["fmod", "remainder"], ["fabs", "absolute"], ["abs", "absolute"],
            ["fmax", "maximum"], ["fmin", "minimum"],
        ]
        for destination, source in equivalences:
            tests.extend([
                [destination, inputs, output] for name, inputs, output in tests 
                if name == source])
        
        for name, inputs, output in tests:
            if name not in ufuncs:
                continue
            self.assertEqual(getattr(numpy, name)(*inputs), output)
        
        untested = [x for x in ufuncs if x not in [t[0] for t in tests]]
        logging.warning(
            "The following ufuncs were not tested: {}".format(
                ", ".join(untested)))

if __name__ == "__main__":
    unittest.main()

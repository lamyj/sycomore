import logging
import math
import pickle
import unittest

import numpy

import sycomore

from test_case import TestCase

class TestArrayQ(TestCase):
    def test_creation(self):
        q = sycomore.ArrayQ([
            [1*sycomore.units.m, 2*sycomore.units.m],
            [3*sycomore.units.m, 4*sycomore.units.m]])
        self.assertEqual(q, sycomore.ArrayQ([[1, 2], [3, 4]], sycomore.Length))
        
        q = sycomore.ArrayQ(numpy.array([
            [1*sycomore.units.m, 2*sycomore.units.m],
            [3*sycomore.units.m, 4*sycomore.units.m]]))
        self.assertEqual(q, sycomore.ArrayQ([[1, 2], [3, 4]], sycomore.Length))
        
        q = sycomore.units.m * [[1, 2], [3, 4]]
        self.assertEqual(q, sycomore.ArrayQ([[1, 2], [3, 4]], sycomore.Length))
        
        q = [[1, 2], [3, 4]] * sycomore.units.m
        self.assertEqual(q, sycomore.ArrayQ([[1, 2], [3, 4]], sycomore.Length))
        
        q = sycomore.ArrayQ(numpy.array([[1, 2], [3, 4]]) * sycomore.units.m)
        self.assertEqual(q, sycomore.ArrayQ([[1, 2], [3, 4]], sycomore.Length))
    
    def test_len(self):
        q = sycomore.ArrayQ([
            [1*sycomore.units.m, 2*sycomore.units.m, 3*sycomore.units.m],
            [4*sycomore.units.m, 5*sycomore.units.m, 6*sycomore.units.m]])
        self.assertEqual(len(q), 2)
    
    def test_shape(self):
        q = sycomore.ArrayQ([
            [1*sycomore.units.m, 2*sycomore.units.m, 3*sycomore.units.m],
            [4*sycomore.units.m, 5*sycomore.units.m, 6*sycomore.units.m]])
        self.assertEqual(q.shape, (2, 3))
    
    def test_comparison(self):
        q1 = [[1, 2], [3, 4]] * sycomore.units.m
        q2 = [[0, 2], [3, 4]] * sycomore.units.m
        q3 = [[1, 2], [3, 4]] * sycomore.units.s
        scalar = sycomore.ArrayQ([[1, 2], [3, 4]])
        
        self.assertTrue(q1 == q1)
        self.assertFalse(q1 == q2)
        self.assertFalse(q1 == q3)
        self.assertTrue(scalar == [[1, 2], [3, 4]])
        self.assertTrue([[1, 2], [3, 4]] == scalar)
        self.assertFalse(scalar == [[0, 2], [3, 4]])
        self.assertFalse([[0, 2], [3, 4]] == scalar)

        self.assertFalse(q1 != q1)
        self.assertTrue(q1 != q2)
        self.assertTrue(q1 != q3)
        self.assertFalse(scalar != [[1, 2], [3, 4]])
        self.assertFalse([[1, 2], [3, 4]] != scalar)
        self.assertTrue(scalar != [[0, 2], [3, 4]])
        self.assertTrue([[0, 2], [3, 4]] !=scalar)
    
    def test_getitem(self):
        from sycomore.units import m
        
        q = [[1, 2], [3, 4]] * m
        self.assertEqual(q[1, 0], 3 * m)
        self.assertEqual(q[-1, -2], 3 * m)
        with self.assertRaises(IndexError):
            q[1]
        with self.assertRaises(IndexError):
            q[2, 0]
        
        q = sycomore.Vector3Q([1*m, 2*m, 3*m])
        self.assertEqual(q[1], 2 * m)
        self.assertEqual(q[-2], 2 * m)
        with self.assertRaises(IndexError):
            q[0, 1]
        with self.assertRaises(IndexError):
            q[3]
    
    def test_setitem(self):
        from sycomore.units import m
        
        q = [[1, 2], [3, 4]] * m
        
        q[1, 0] = 42 * m
        self.assertEqual(q[1, 0], 42 * m)
        
        q[-1, -2] = 43 * m
        self.assertEqual(q[1, 0], 43 * m)
        
        with self.assertRaises(IndexError):
            q[1] = 1 * m
        with self.assertRaises(IndexError):
            q[2, 0] = 1 * m
        
        q = sycomore.Vector3Q([1*m, 2*m, 3*m])
        
        q[1] = 42 * m
        self.assertEqual(q[1], 42 * sycomore.units.m)
        
        q[-2] = 43 * m
        self.assertEqual(q[1], 43 * sycomore.units.m)
        with self.assertRaises(IndexError):
            q[0, 1] = 1 * m
        with self.assertRaises(IndexError):
            q[3] = 1 * m
        
    def test_addition_in_place(self):
        q = [[1, 2], [3, 4]] * sycomore.units.m
        q += [[5, 6], [7, 8]] * sycomore.units.m
        self.assertEqual(q, [[6, 8], [10, 12]] * sycomore.units.m)
        
        q += 9 * sycomore.units.m
        self.assertEqual(q, [[15, 17], [19, 21]] * sycomore.units.m)
        
        with self.assertRaises(RuntimeError):
            q += 2 * sycomore.units.s
        with self.assertRaises(RuntimeError):
            q += 3
        
        q = sycomore.ArrayQ([[1, 2], [3, 4]])
        q += [[5, 6], [7, 8]]
        self.assertEqual(q, sycomore.ArrayQ([[6, 8], [10, 12]]))
        
        q += 9
        self.assertEqual(q, sycomore.ArrayQ([[15, 17], [19, 21]]))
    
    def test_subtraction_in_place(self):
        q = [[1, 2], [3, 4]] * sycomore.units.m
        q -= [[8, 7], [6, 5]] * sycomore.units.m
        self.assertEqual(q, [[-7, -5], [-3, -1]] * sycomore.units.m)
        
        q -= 9 * sycomore.units.m
        self.assertEqual(q, [[-16, -14], [-12, -10]] * sycomore.units.m)
        
        with self.assertRaises(RuntimeError):
            q -= 2 * sycomore.units.s
        with self.assertRaises(RuntimeError):
            q -= 3
        
        q = sycomore.ArrayQ([[1, 2], [3, 4]])
        q -= [[8, 7], [6, 5]]
        self.assertEqual(q, sycomore.ArrayQ([[-7, -5], [-3, -1]]))
        
        q -= 9
        self.assertEqual(q, sycomore.ArrayQ([[-16, -14], [-12, -10]]))
    
    def test_multiplication_in_place(self):
        q = [[1, 2], [3, 4]] * sycomore.units.m
        q *= [[5, 6], [7, 8]] * sycomore.units.s
        self.assertEqual(
            q, [[5, 12], [21, 32]] * sycomore.units.m*sycomore.units.s)
        
        q *= 2 * sycomore.units.kg
        self.assertEqual(
            q, [[10, 24], [42, 64]] * sycomore.units.m*sycomore.units.s*sycomore.units.kg)
        
        q *= 0.5
        self.assertEqual(
            q, [[5, 12], [21, 32]] * sycomore.units.m*sycomore.units.s*sycomore.units.kg)
    
    def test_division_in_place(self):
        q = [[10, 24], [42, 64]] * sycomore.units.m*sycomore.units.s*sycomore.units.kg
        q /= [[5, 6], [7, 8]] * sycomore.units.kg
        self.assertEqual(
            q, [[2, 4], [6, 8]] * sycomore.units.m*sycomore.units.s)
        
        q /= 2 * sycomore.units.s
        self.assertEqual(q, [[1, 2], [3, 4]] * sycomore.units.m)
        
        q /= 0.5
        self.assertEqual(q, [[2, 4], [6, 8]] * sycomore.units.m)
    
    def test_floordiv_in_place(self):
        q = [[7, 11], [13, 17]] * sycomore.units.m*sycomore.units.s
        q //= [[5, 4], [3, 2]] * sycomore.units.s
        self.assertEqual(q, [[1, 2], [4, 8]] * sycomore.units.m)
        
        q = [[7, 11], [13, 17]] * sycomore.units.m*sycomore.units.s
        q //= 2 * sycomore.units.s
        self.assertEqual(q, [[3, 5], [6, 8]] * sycomore.units.m)
        
        q //= 2
        self.assertEqual(q, [[1, 2], [3, 4]] * sycomore.units.m)
    
    def test_modulo_in_place(self):
        q = [[7, 11], [13, 17]] * sycomore.units.m
        q %= [[5, 4], [3, 2]] * sycomore.units.m
        self.assertEqual(q, [[2, 3], [1, 1]] * sycomore.units.m)
        
        q = [[7, 10], [5, 8]] * sycomore.units.m
        q %= 4 * sycomore.units.m
        self.assertEqual(q, [[3, 2], [1, 0]] * sycomore.units.m)
        
        q = [[7, 10], [5, 8]] * sycomore.units.m
        q %= 4
        self.assertEqual(q, [[3, 2], [1, 0]] * sycomore.units.m)
        
        with self.assertRaises(RuntimeError):
            q %= [[5, 4], [3, 2]] * sycomore.units.s
        with self.assertRaises(RuntimeError):
            q %= 2 * sycomore.units.s
    
    def test_convert_to(self):
        q = [[1, 2], [3, 4]] * sycomore.units.cm
        numpy.testing.assert_allclose(
            q.convert_to(sycomore.units.mm), [[10, 20], [30, 40]])
        with self.assertRaises(RuntimeError):
            q.convert_to(sycomore.units.s)
    
    def test_unary_plus(self):
        q = [[1, 2], [3, 4]] * sycomore.units.cm
        self.assertTrue(+q == q)
    
    def test_unary_minus(self):
        q = [[1, 2], [3, 4]] * sycomore.units.cm
        self.assertTrue(-q == -1*q)
    
    def test_addition(self):
        q1 = [[1, 2], [3, 4]] * sycomore.units.m
        q2 = [[5, 6], [7, 8]] * sycomore.units.m
        self.assertEqual(q1 + q2, [[6, 8], [10, 12]] * sycomore.units.m)
        
        q2 = 9 * sycomore.units.m
        self.assertEqual(q1 + q2, [[10, 11], [12, 13]] * sycomore.units.m)
        self.assertEqual(q2 + q1, [[10, 11], [12, 13]] * sycomore.units.m)
        
        with self.assertRaises(RuntimeError):
            q1 + 2 * sycomore.units.s
        with self.assertRaises(RuntimeError):
            2 * sycomore.units.s + q1
        with self.assertRaises(RuntimeError):
            q1 + 3
        with self.assertRaises(RuntimeError):
            3 + q1
        
        q1 = sycomore.ArrayQ([[1, 2], [3, 4]])
        q2 = [[5, 6], [7, 8]]
        self.assertEqual(q1 + q2, sycomore.ArrayQ([[6, 8], [10, 12]]))
        self.assertEqual(q2 + q1, sycomore.ArrayQ([[6, 8], [10, 12]]))
        
        q2 = 9
        self.assertEqual(q1 + q2, sycomore.ArrayQ([[10, 11], [12, 13]]))
        self.assertEqual(q2 + q1, sycomore.ArrayQ([[10, 11], [12, 13]]))
    
    def test_subtraction(self):
        q1 = [[1, 2], [3, 4]] * sycomore.units.m
        q2 = [[8, 7], [6, 5]] * sycomore.units.m
        self.assertEqual(q1 - q2, [[-7, -5], [-3, -1]] * sycomore.units.m)
        
        q2 = 9 * sycomore.units.m
        self.assertEqual(q1 - q2, [[-8, -7], [-6, -5]] * sycomore.units.m)
        self.assertEqual(-(q2 - q1), [[-8, -7], [-6, -5]] * sycomore.units.m)
        
        with self.assertRaises(RuntimeError):
            q1 - 2 * sycomore.units.s
        with self.assertRaises(RuntimeError):
            2 * sycomore.units.s - q1
        with self.assertRaises(RuntimeError):
            q1 - 3
        with self.assertRaises(RuntimeError):
            3 - q1
        
        q1 = sycomore.ArrayQ([[1, 2], [3, 4]])
        q2 = [[8, 7], [6, 5]]
        self.assertEqual(q1 - q2, sycomore.ArrayQ([[-7, -5], [-3, -1]]))
        self.assertEqual(-(q2 - q1), sycomore.ArrayQ([[-7, -5], [-3, -1]]))
        
        q2 = 9
        self.assertEqual(q1 - q2, [[-8, -7], [-6, -5]])
        self.assertEqual(-(q2 - q1), [[-8, -7], [-6, -5]])
    
    def test_multiplication(self):
        q1 = [[1, 2], [3, 4]] * sycomore.units.m
        q2 = [[5, 6], [7, 8]] * sycomore.units.s
        self.assertEqual(
            q1 * q2, [[5, 12], [21, 32]] * sycomore.units.m*sycomore.units.s)
        
        q2 = 2 * sycomore.units.s
        self.assertEqual(
            q1 * q2, [[2, 4], [6, 8]] * sycomore.units.m*sycomore.units.s)
        self.assertEqual(
            q2 * q1, [[2, 4], [6, 8]] * sycomore.units.m*sycomore.units.s)
        
        q2 = 2
        self.assertEqual(q1 * q2, [[2, 4], [6, 8]] * sycomore.units.m)
        self.assertEqual(q2 * q1, [[2, 4], [6, 8]] * sycomore.units.m)
    
    def test_division(self):
        q1 = [[10, 24], [42, 64]] * sycomore.units.m*sycomore.units.s
        q2 = [[5, 6], [7, 8]] * sycomore.units.s
        self.assertEqual(q1 / q2, [[2, 4], [6, 8]] * sycomore.units.m)
        
        q2 = 2 * sycomore.units.s
        self.assertEqual(
            q1 / q2, [[5, 12], [21, 32]] * sycomore.units.m)
        self.assertEqual(
            1/(q2 / q1), [[5, 12], [21, 32]] * sycomore.units.m)
        
        q2 = 2
        self.assertEqual(
            q1 / q2, [[5, 12], [21, 32]] * sycomore.units.m*sycomore.units.s)
        self.assertEqual(
            1/(q2 / q1), [[5, 12], [21, 32]] * sycomore.units.m*sycomore.units.s)
    
    def test_floordiv(self):
        q1 = [[7, 11], [13, 17]] * sycomore.units.m*sycomore.units.s
        q2 = [[5, 4], [3, 2]] * sycomore.units.s
        self.assertEqual(q1 // q2, [[1, 2], [4, 8]] * sycomore.units.m)
        
        q1 = [[7, 11], [13, 17]] * sycomore.units.m*sycomore.units.s
        q2 = 2 * sycomore.units.s
        self.assertEqual(q1 // q2, [[3, 5], [6, 8]] * sycomore.units.m)
        
        q1 = [[7, 11], [13, 17]] * sycomore.units.m
        q2 = 2
        self.assertEqual(q1 // q2, [[3, 5], [6, 8]] * sycomore.units.m)
        
        q1 = 17 * sycomore.units.m*sycomore.units.s
        q2 = [[2, 3], [4, 5]] * sycomore.units.s
        self.assertEqual(q1 // q2, [[8, 5], [4, 3]] * sycomore.units.m)
        
        q1 = 17
        q2 = [[2, 3], [4, 5]] * sycomore.units.m
        self.assertEqual(
            q1 // q2, [[8, 5], [4, 3]] / sycomore.units.m)
    
    def test_modulo(self):
        q1 = [[7, 11], [13, 17]] * sycomore.units.m
        q2 = [[5, 4], [3, 2]] * sycomore.units.m
        self.assertEqual(q1 % q2, [[2, 3], [1, 1]] * sycomore.units.m)
        
        q1 = [[7, 10], [5, 8]] * sycomore.units.m
        q2 = 4 * sycomore.units.m
        self.assertEqual(q1 % q2, [[3, 2], [1, 0]] * sycomore.units.m)
        
        q2 = 4
        self.assertEqual(q1 % q2, [[3, 2], [1, 0]] * sycomore.units.m)
        
        with self.assertRaises(RuntimeError):
            q1 % ([[5, 4], [3, 2]] * sycomore.units.s)
        with self.assertRaises(RuntimeError):
            q1 % (2 * sycomore.units.s)
    
    def test_divmod(self):
        q1 = [[7, 11], [13, 17]] * sycomore.units.m
        q2 = [[5, 4], [3, 2]] * sycomore.units.m
        self.assertEqual(divmod(q1, q2), (q1 // q2, q1 % q2))
        
        q2 = [[5, 4], [3, 2]]
        self.assertEqual(divmod(q1, q2), (q1 // q2, q1 % q2))
    
    def test_abs(self):
        self.assertEqual(
            abs([[1, -2], [3, -4]] / sycomore.units.m),
            [[1, 2], [3, 4]] / sycomore.units.m)
    
    def test_pow(self):
        q = [[1, -2], [3, -4]] / sycomore.units.m
        self.assertEqual(q**3, [[1, -8], [27, -64]] / sycomore.units.m**3)
        self.assertEqual(pow(q, 3), [[1, -8], [27, -64]] / sycomore.units.m**3)
    
    def test_round(self):
        q = [[9.2, -9.7], [-9, 0]] * pow(sycomore.units.m, 1.5)
        self.assertEqual(
            round(q), [[9, -10], [-9, 0]] * pow(sycomore.units.m, 1.5))
    
    def test_trunc(self):
        q = [[9.2, -9.7], [-9, 0]] * pow(sycomore.units.m, 1.5)
        self.assertEqual(
            math.trunc(q), [[9, -9], [-9, 0]] * pow(sycomore.units.m, 1.5))
    
    def test_floor(self):
        q = [[9.2, -9.7], [-9, 0]] * pow(sycomore.units.m, 1.5)
        self.assertEqual(
            math.floor(q), [[9, -10], [-9, 0]] * pow(sycomore.units.m, 1.5))
    
    def test_ceil(self):
        q = [[9.2, -9.7], [-9, 0]] * pow(sycomore.units.m, 1.5)
        self.assertEqual(
            math.ceil(q), [[10, -9], [-9, 0]] * pow(sycomore.units.m, 1.5))
    
    def test_pickle(self):
        q = sycomore.ArrayQ([-1, 0.5], sycomore.Dimensions(7,6,5,4,3,2,1))
        self.assertEqual(pickle.loads(pickle.dumps(q)), q)
    
    def test_hash(self):
        quantities = set()

        quantities.add([[1, 2], [3, 4]] * sycomore.units.m)
        quantities.add([[1, 2], [3, 4]] * sycomore.units.s)
        quantities.add([[0, 2], [3, 4]] * sycomore.units.m)
        self.assertEqual(len(quantities), 3)
        
        quantities.add([[1, 2], [3, 4]] * sycomore.units.m)
        self.assertEqual(len(quantities), 3)
        
    def test_ufuncs(self):
        ufuncs = [
            x for x in dir(numpy) if isinstance(getattr(numpy, x), numpy.ufunc)]
        non_object_ufuncs = [
            x for x in ufuncs if not any(
                t.endswith("O->O") for t in getattr(numpy, x).types)
        ]
        logging.info(
            "The following ufuncs do not operate on objects: {}".format(
                ", ".join(non_object_ufuncs)))
        ufuncs = [x for x in ufuncs if x not in non_object_ufuncs]
        
        not_applicable_ufuncs = [
            "greater", "greater_equal", "less", "less_equal",
            "equal", "not_equal", "maximum", "minimum", "fmax", "fmin", "sign",
            
            "conj", "conjugate", # Quantities are real-valued
            "deg2rad", "degrees", "rad2deg", "radians", 
            "bitwise_and", "bitwise_or", "bitwise_xor", "bitwise_not",
            "bitwise_count", "bitwise_invert", "bitwise_left_shift",
            "bitwise_right_shift",
            "invert", "left_shift", "right_shift", "logical_and", "logical_or", 
            "logical_xor", "logical_not", "isnat", 
            "gcd", "lcm", "matmul", "matvec", "vecdot", "vecmat"
        ]
        ufuncs = [x for x in ufuncs if x not in not_applicable_ufuncs]
        
        from sycomore.units import m, deg
        Scalar = sycomore.Quantity(1, sycomore.Dimensions())
        tests = [
            ["add", [[1, 2]*m, [2]*m], [3, 4]*m],
            ["add", [[1, 2]*Scalar, 2], [3, 4]*Scalar],
            ["add", [1, [2, 3]*Scalar], [3, 4]*Scalar],
            
            ["subtract", [[3, 4]*m, 2*m], [1, 2]*m],
            ["subtract", [[3, 4]*Scalar, 2], [1, 2]*Scalar],
            ["subtract", [2, [3, 4]*Scalar], [-1, -2]*Scalar],
            
            ["multiply", [[2, 4]*m, 3*m], [6, 12]*m**2],
            ["multiply", [[2, 4]*m, 3], [6, 12]*m],
            ["multiply", [3, [2, 4]*m], [6, 12]*m],
            
            ["divide", [[3, 2]*m, 2*m], [1.5, 1]*Scalar],
            ["divide", [[3, 2]*m, 2], [1.5, 1]*m],
            ["divide", [3, [1, 2]*m], [3, 1.5]/m],
            # Same tests for true_divide
            
            ["floor_divide", [[4, 3]*m, 2*m], [2, 1]*Scalar],
            ["floor_divide", [[4, 3]*m, 2], [2, 1]*m],
            ["floor_divide", [4, [3, 2]*m], [1, 2]/m],
            
            ["negative", [[-1, +1]*m], [+1, -1]*m],
            ["positive", [[-1, +1]*m], [-1, +1]*m],
            ["power", [[2, 3]*m, 2], [4, 9]*m**2],
            
            ["remainder", [[3, 2]*m, 2], [1, 0]*m],
            ["remainder", [[3, 2]*m, 2*m], [1, 0]*m],
            # Same tests for mod and fmod
            
            ["absolute", [[-1, 1]*m], [1, 1]*m],
            # Same tests for fabs and abs
            
            ["rint", [[-1.7, -2.2]*m], [-2, -2]*m],
                        
            ["exp", [[0, 1]*Scalar], [1, 2.718281828]*Scalar],
            ["exp2", [[0, 1]*Scalar], [1, 2]*Scalar],
            ["log", [[1, 2.718281828]*Scalar], [0, 1]*Scalar],
            ["log2", [[4, 8]*Scalar], [2, 3]*Scalar],
            ["log10", [[100, 1000]*Scalar], [2, 3]*Scalar],
            ["expm1", [[0, 1]*Scalar], [0, 1.718281828]*Scalar],
            ["log1p", [[0, 1.718281828]*Scalar], [0, 1]*Scalar],
            
            ["sqrt", [[4, 9]*m**2], [2, 3]*m],
            ["square", [[2, 3]*m], [4, 9]*m**2],
            ["cbrt", [[8, 27]*m**3], [2, 3]*m],
            ["reciprocal", [[2, 4]*m], [0.5, 0.25]/m],
            
            ["sin", [[0, 90]*deg], [0, 1]*Scalar],
            ["cos", [[0, 90]*deg], [1, 0]*Scalar],
            ["tan", [[0, 45]*deg], [0, 1]*Scalar],
            ["arcsin", [[0, 0.5]*Scalar], [0, 30]*deg],
            ["arccos", [[1, 0.5]*Scalar], [0, 60]*deg],
            ["arctan", [[0, 1]*Scalar], [0, 45]*deg],
            ["arctan2", [[0, 1, 2]*Scalar, [1, 1, 2]*Scalar], [0, 45, 45]*deg],
            ["hypot", [[3, 5]*Scalar, [4, 12]*Scalar], [5, 13]*Scalar],
            
            ["sinh", [[0, 1]*Scalar], [0, numpy.sinh(1)]*Scalar],
            ["cosh", [[0, 1]*Scalar], [1, numpy.cosh(1)]*Scalar],
            ["tanh", [[0, 1]*Scalar], [0, numpy.tanh(1)]*Scalar],
            ["arcsinh", [[0, 1]*Scalar], [0, numpy.arcsinh(1)]*Scalar],
            ["arccosh", [[1, 2]*Scalar], [0, numpy.arccosh(2)]*Scalar],
            ["arctanh", [[0, 0.5]*Scalar], [0, numpy.arctanh(0.5)]*Scalar],
            
            ["ceil", [[-2.1, -1.9]*m], [-2, -1]*m],
            ["floor", [[-2.1, -1.9]*m], [-3, -2]*m],
            ["trunc", [[-2.1, -1.9]*m], [-2, -1]*m],
        ]
        equivalences = [
            ["true_divide", "divide"], ["mod", "remainder"], 
            ["fmod", "remainder"], ["fabs", "absolute"], ["abs", "absolute"],
            ["fmax", "maximum"], ["fmin", "minimum"], ["pow", "power"],
            *([f"a{name}", f"arc{name}"] for name in [
                "cos", "sin", "tan", "tan2", "sinh", "cosh", "tanh"])
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
        if untested:
            logging.warning(
                "The following ufuncs were not tested: {}".format(
                    ", ".join(untested)))
    
if __name__ == "__main__":
    unittest.main()

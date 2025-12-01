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
    
    def test_iter(self):
        qa = sycomore.ArrayQ([
            [1*sycomore.units.m, 2*sycomore.units.m, 3*sycomore.units.m],
            [4*sycomore.units.m, 5*sycomore.units.m, 6*sycomore.units.m]])
        
        ql = [q for q in qa]
        self.assertEqual(ql, [i * sycomore.units.m for i in range(1, 7)])
        
        ql = list(qa)
        self.assertEqual(ql, [i * sycomore.units.m for i in range(1, 7)])
    
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
        with self.assertRaises(RuntimeError):
            q.scalar
        
        q = sycomore.ArrayQ([[1, 2], [3, 4]])
        numpy.testing.assert_allclose(q.scalar, q.magnitude)
    
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
    
if __name__ == "__main__":
    unittest.main()

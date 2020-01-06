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

        self.assertTrue(q1 == q1)
        self.assertFalse(q1 == q2)
        self.assertFalse(q1 == q3)

        self.assertFalse(q1 != q1)
        self.assertTrue(q1 != q2)
        self.assertTrue(q1 != q3)

        self.assertTrue(q1 < q2)
        self.assertFalse(q2 <= q1)
        self.assertFalse(q1 > q2)
        self.assertFalse(q1 >= q2)
        with self.assertRaises(Exception):
            q1 < q3
        with self.assertRaises(Exception):
            q1 <= q3
        with self.assertRaises(Exception):
            q1 > q3
        with self.assertRaises(Exception):
            q1 >= q3

    def test_addition_in_place(self):
        q1 = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(5, sycomore.Dimensions(1,0,0,0,0,0,0))
        q1 += q2
        self.assertEqual(q1, r)

        q3 = sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0))
        with self.assertRaises(Exception):
            q1 += q3

    def test_subtraction_in_place(self):
        q1 = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(-1, sycomore.Dimensions(1,0,0,0,0,0,0))
        q1 -= q2
        self.assertEqual(q1, r)

        q3 = sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0))
        with self.assertRaises(Exception):
            q1 -= q3

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
        r = sycomore.Quantity(5, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(q1+q2, r)

        q3 = sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0))
        with self.assertRaises(Exception):
            q1 + q3

    def test_subtraction(self):
        q1 = sycomore.Quantity(2, sycomore.Dimensions(1,0,0,0,0,0,0))
        q2 = sycomore.Quantity(3, sycomore.Dimensions(1,0,0,0,0,0,0))
        r = sycomore.Quantity(-1, sycomore.Dimensions(1,0,0,0,0,0,0))
        self.assertEqual(q1-q2, r)

        q3 = sycomore.Quantity(2, sycomore.Dimensions(0,1,0,0,0,0,0))
        with self.assertRaises(Exception):
            q1 - q3

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
        quantity_functions = [
            "add", "subtract", "multiply", "divide", 
            # TODO? logaddexp, logaddexp2 
            "true_divide", "floor_divide", "negative", "positive",
            "remainder", "mod", "fmod", 
            # TODO divmod
            "absolute", "fabs",
            # TODO? sign, heaviside
            # Not applicable: conj, conjugate
            "sqrt", "square", "cbrt", "reciprocal",
            # Not applicable: deg2rad, rad2deg
            # Not applicable: bitwise_and, bitwise_or, bitwise_xor, invert, 
            # left_shift, right_shift
            # Not applicable: logical_and, logical_or, logical_xor, logical_not
            "maximum", "minimum", "fmax", "fmin",
            # TODO? isfinite, isinf
            # Not applicable: isnat
            # TODO? signbit, copysign, nextafter, spacing, modf, ldexp, frexp
            "floor", "ceil", "trunc"
        ]
        quantity_functions = list(zip(
            quantity_functions, len(quantity_functions)*[False]))
        
        scalar_only_functions = [
            "exp", "exp2", "log", "log2", "log10", "expm1", "log1p",
            # TODO? "gcd", "lcm",
            "sin", "cos", "tan", "arcsin", "arccos", "arctan", "arctan2",
        ]
        scalar_only_functions = list(zip(
            scalar_only_functions, len(scalar_only_functions)*[True]))
        
        Scalar = sycomore.Quantity(1, sycomore.Dimensions())
        quantities = [0.123*Scalar, 3*Scalar]
        for name, scalar_only in quantity_functions+scalar_only_functions:
            function = getattr(numpy, name, None)
            if function is None:
                continue
            arguments = quantities[:function.nin]
            scalar_arguments = [q.magnitude for q in arguments]
            result = function(*arguments)
            self.assertEqual(result, Scalar*function(*scalar_arguments))
            
            non_scalar_arguments = [q*sycomore.units.m for q in arguments]
            if not scalar_only:
                function(*non_scalar_arguments)
            else:
                with self.assertRaises(Exception):
                    function(*non_scalar_arguments)
        
        comparators = [
            "greater", "greater_equal", "less", "less_equal", 
            "not_equal", "equal"]
        for name in comparators:
            function = getattr(numpy, name)
            self.assertEqual(
                function(*quantities[:2]),
                function(*[q.magnitude for q in quantities[:2]]))
        
        # power has a different signature
        result = numpy.power(4*sycomore.units.m**6, 0.5)
        self.assertEqual(result, 2*sycomore.units.m**3)

if __name__ == "__main__":
    unittest.main()

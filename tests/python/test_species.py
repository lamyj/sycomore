import math
import unittest

import sycomore
from sycomore.units import *

class TestSpecies(unittest.TestCase):
    def test_quantity_constructor_full(self):
        species = sycomore.Species(
            1*ms, 0.01*kHz, 3*um*um/ms, 1./(15*ds), 0.9*rad/s, 1.1)
        self.assertEqual(species.R1, 1.*kHz)
        self.assertEqual(species.T1, 1.*ms)
        self.assertEqual(species.R2, 10.*Hz)
        self.assertEqual(species.T2, 0.1*s)
        self.assertEqual(species.D[0], 3e-9*m*m/s)
        self.assertEqual(species.R2_prime, (1/1.5)*Hz)
        self.assertEqual(species.T2_prime, 1.5*s)
        self.assertEqual(species.delta_omega, 0.9*rad/s)
        self.assertEqual(species.w, 1.1)

    def test_quantity_constructor_partial(self):
        species = sycomore.Species(
            1*ms, 0.01*kHz, R2_prime=1./(15*ds))
        self.assertEqual(species.R1, 1.*kHz)
        self.assertEqual(species.T1, 1.*ms)
        self.assertEqual(species.R2, 10.*Hz)
        self.assertEqual(species.T2, 0.1*s)
        self.assertEqual(species.D[0], 0*m*m/s)
        self.assertEqual(species.R2_prime, (1/1.5)*Hz)
        self.assertEqual(species.T2_prime, 1.5*s)
        self.assertEqual(species.delta_omega, 0*rad/s)
        self.assertEqual(species.w, 1)

    def test_D_scalar(self):
        D_scalar = 1*um*um/ms
        D_tensor = [
            1*um*um/ms, 0*um*um/ms, 0*um*um/ms,
            0*um*um/ms, 1*um*um/ms, 0*um*um/ms,
            0*um*um/ms, 0*um*um/ms, 1*um*um/ms]

        species = sycomore.Species(1*ms, 100*ms, D_scalar)
        self.assertSequenceEqual(species.D, D_tensor)

        species = sycomore.Species(1*ms, 100*ms)
        species.D = D_scalar
        self.assertSequenceEqual(species.D, D_tensor)

    def test_D_tensor(self):

        D = [
            1*um*um/ms, 4*um*um/ms, 7*um*um/ms,
            2*um*um/ms, 5*um*um/ms, 8*um*um/ms,
            3*um*um/ms, 6*um*um/ms, 9*um*um/ms]

        species = sycomore.Species(
            1*ms, 100*ms, sycomore.Array[sycomore.Quantity](D))
        self.assertSequenceEqual(species.D, D)

        species = sycomore.Species(1*ms, 100*ms, D)
        self.assertSequenceEqual(species.D, D)

        species = sycomore.Species(1*ms, 100*ms)
        species.D = sycomore.Array[sycomore.Quantity](D)
        self.assertSequenceEqual(species.D, D)

        species.D = D
        self.assertSequenceEqual(species.D, D)

if __name__ == "__main__":
    unittest.main()

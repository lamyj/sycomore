import math
import unittest

import sycomore
from sycomore.units import *

class TestSpecies(unittest.TestCase):
    def test_scalar_constructor_full(self):
        species = sycomore.Species(1., 10., 3e-9, 1/1.5, 0.9, 1.1)
        self.assertEqual(species.R1, 1.)
        self.assertEqual(species.R2, 10.)
        self.assertEqual(species.D, 3e-9)
        self.assertEqual(species.R2_prime, 1/1.5)
        self.assertEqual(species.delta_omega, 0.9)
        self.assertEqual(species.w, 1.1)

    def test_scalar_constructor_partial(self):
        species = sycomore.Species(1., 10., R2_prime=1/1.5)
        self.assertEqual(species.R1, 1.)
        self.assertEqual(species.R2, 10.)
        self.assertEqual(species.D, 0)
        self.assertEqual(species.R2_prime, 1/1.5)
        self.assertEqual(species.delta_omega, 0.)
        self.assertEqual(species.w, 1.)

    def test_quantity_constructor_full(self):
        species = sycomore.Species(
            1000*ms, 0.01*kHz, 3*um*um/ms, 1./(15*ds), 0.9*rad/s, 1.1)
        self.assertEqual(species.R1, 1.)
        self.assertEqual(species.R2, 10.)
        self.assertEqual(species.D, 3e-9)
        self.assertEqual(species.R2_prime, 1/1.5)
        self.assertEqual(species.delta_omega, 0.9)
        self.assertEqual(species.w, 1.1)

    def test_quantity_constructor_partial(self):
        species = sycomore.Species(
            1000*ms, 0.01*kHz, R2_prime=1./(15*ds))
        self.assertEqual(species.R1, 1.)
        self.assertEqual(species.R2, 10.)
        self.assertEqual(species.D, 0)
        self.assertEqual(species.R2_prime, 1/1.5)
        self.assertEqual(species.delta_omega, 0)
        self.assertEqual(species.w, 1)

if __name__ == "__main__":
    unittest.main()

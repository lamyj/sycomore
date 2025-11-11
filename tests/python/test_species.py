import math
import pickle
import unittest

import numpy

import sycomore
from sycomore.units import *

from test_case import TestCase

class TestSpecies(TestCase):
    def test_quantity_constructor_full(self):
        species = sycomore.Species(1*ms, 0.01*kHz, 3*um*um/ms, 0.9*rad/s)
        self.assertEqual(species.R1, 1.*kHz)
        self.assertEqual(species.T1, 1.*ms)
        self.assertEqual(species.R2, 10.*Hz)
        self.assertEqual(species.T2, 0.1*s)
        self.assertEqual(species.D[0,0], 3e-9*m*m/s)
        self.assertEqual(species.delta_omega, 0.9*rad/s)
    
    def test_quantity_constructor_partial(self):
        species = sycomore.Species(1*ms, 0.01*kHz, delta_omega=0.9*rad/s)
        self.assertEqual(species.R1, 1.*kHz)
        self.assertEqual(species.T1, 1.*ms)
        self.assertEqual(species.R2, 10.*Hz)
        self.assertEqual(species.T2, 0.1*s)
        self.assertEqual(species.D[0,0], 0*m*m/s)
        self.assertEqual(species.delta_omega, 0.9*rad/s)
    
    def test_D_scalar(self):
        D_scalar = 1*um*um/ms
        D_tensor = sycomore.Matrix3x3Q([
            [1*um*um/ms, 0*um*um/ms, 0*um*um/ms],
            [0*um*um/ms, 1*um*um/ms, 0*um*um/ms],
            [0*um*um/ms, 0*um*um/ms, 1*um*um/ms]])
        
        species = sycomore.Species(1*ms, 100*ms, D_scalar)
        self.assertEqual(species.D, D_tensor)

        species = sycomore.Species(1*ms, 100*ms)
        species.D = D_scalar
        self.assertEqual(species.D, D_tensor)

    def test_D_tensor(self):

        D = sycomore.Matrix3x3Q([
            [1*um*um/ms, 4*um*um/ms, 7*um*um/ms],
            [2*um*um/ms, 5*um*um/ms, 8*um*um/ms],
            [3*um*um/ms, 6*um*um/ms, 9*um*um/ms]])

        species = sycomore.Species(1*ms, 100*ms, D)
        self.assertEqual(species.D, D)

        species = sycomore.Species(1*ms, 100*ms, D)
        self.assertEqual(species.D, D)

        species = sycomore.Species(1*ms, 100*ms)
        species.D = D
        self.assertEqual(species.D, D)

        species.D = D
        self.assertEqual(species.D, D)
    
    def test_pickle(self):
        D = sycomore.Matrix3x3Q([
            [1*um*um/ms, 4*um*um/ms, 7*um*um/ms],
            [2*um*um/ms, 5*um*um/ms, 8*um*um/ms],
            [3*um*um/ms, 6*um*um/ms, 9*um*um/ms]])

        species = sycomore.Species(1*ms, 100*ms, D, 0.9*rad/s)
        
        other_species = pickle.loads(pickle.dumps(species))
        self.assertEqual(species.T1, other_species.T1)
        self.assertEqual(species.T2, other_species.T2)
        self.assertEqual(species.D, other_species.D)
        self.assertEqual(species.delta_omega, other_species.delta_omega)

if __name__ == "__main__":
    unittest.main()

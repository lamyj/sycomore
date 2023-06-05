import math
import pickle
import unittest

import numpy

import sycomore
from sycomore.units import *

class TestSpecies(unittest.TestCase):
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
        D_tensor = [
            [1*um*um/ms, 0*um*um/ms, 0*um*um/ms],
            [0*um*um/ms, 1*um*um/ms, 0*um*um/ms],
            [0*um*um/ms, 0*um*um/ms, 1*um*um/ms]]

        species = sycomore.Species(1*ms, 100*ms, D_scalar)
        self._test_quantity_array(species.D, D_tensor)

        species = sycomore.Species(1*ms, 100*ms)
        species.D = D_scalar
        self._test_quantity_array(species.D, D_tensor)

    def test_D_tensor(self):

        D = [
            [1*um*um/ms, 4*um*um/ms, 7*um*um/ms],
            [2*um*um/ms, 5*um*um/ms, 8*um*um/ms],
            [3*um*um/ms, 6*um*um/ms, 9*um*um/ms]]

        species = sycomore.Species(1*ms, 100*ms, D)
        self._test_quantity_array(species.D, D)

        species = sycomore.Species(1*ms, 100*ms, D)
        self._test_quantity_array(species.D, D)

        species = sycomore.Species(1*ms, 100*ms)
        species.D = numpy.array(D)
        self._test_quantity_array(species.D, D)

        species.D = D
        self._test_quantity_array(species.D, D)
    
    def test_pickle(self):
        D = [
            [1*um*um/ms, 4*um*um/ms, 7*um*um/ms],
            [2*um*um/ms, 5*um*um/ms, 8*um*um/ms],
            [3*um*um/ms, 6*um*um/ms, 9*um*um/ms]]

        species = sycomore.Species(1*ms, 100*ms, D, 0.9*rad/s)
        
        other_species = pickle.loads(pickle.dumps(species))
        self.assertEqual(species.T1, other_species.T1)
        self.assertEqual(species.T2, other_species.T2)
        self._test_quantity_array(species.D, other_species.D)
        self.assertEqual(species.delta_omega, other_species.delta_omega)
    
    def _test_quantity_array(self, left, right):
        self.assertEqual(numpy.shape(left), numpy.shape(right))
        self.assertSequenceEqual(
            [x.dimensions for x in numpy.ravel(left)],
            [x.dimensions for x in numpy.ravel(right)])
        self.assertSequenceEqual(
            [x.magnitude for x in numpy.ravel(left)],
            [x.magnitude for x in numpy.ravel(right)])

if __name__ == "__main__":
    unittest.main()
